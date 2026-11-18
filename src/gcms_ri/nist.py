"""NIST MS Search background integration and report parsing."""

import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time


_HIT = re.compile(
    r"Hit\s+(\d+)\s*:\s*<<(.*?)>>;\s*<<(.*?)>>;\s*MF:\s*(\d+);"
    r"\s*RMF:\s*(\d+);\s*Prob:\s*([\d.]+);\s*CAS:\s*(.*?);"
    r"\s*Mw:\s*(.*?);\s*Lib:\s*<<(.*?)>>;\s*Id:\s*(.*?);\s*RI:\s*([^\r\n]*?)(?:\.)?\s*$"
)


def parse_report(path):
    data = path.read_bytes()
    try:
        content = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        content = data.decode("mbcs" if os.name == "nt" else "cp1252", errors="replace")
    results = {}
    for block in re.split(r"(?m)^Unknown:\s+", content)[1:]:
        query = block.split()[0]
        if query in results:
            raise ValueError(f"Duplicate NIST report query: {query}")
        hits = []
        for line in block.splitlines()[1:]:
            if not re.match(r"^\s*Hit\s+\d", line):
                continue
            match = _HIT.match(line.lstrip())
            if not match:
                raise ValueError(f"Unrecognized NIST hit format for {query}: {line}")
            rank, name, formula, mf, rmf, prob, cas, mw, library, number, ri = match.groups()
            hits.append(dict(rank=int(rank), name=name, formula=formula, mf=int(mf), rmf=int(rmf),
                             prob=float(prob), cas=cas.strip(), mw=mw.strip(), library=library,
                             library_id=number.strip(), library_ri=ri.strip()))
        results[query] = hits
    if not results:
        raise ValueError("NIST report contains no Unknown: query blocks.")
    return results


def nist_ready(directory):
    executable = directory / "nistms$.exe"
    if os.name != "nt" or not executable.is_file() or not (directory / "nistms.INI").is_file():
        raise ValueError("Supply a Windows NIST MS Search directory containing nistms$.exe and nistms.INI.")
    return executable


def search(msp, directory, output, timeout=600, temp_root=None):
    executable = nist_ready(directory)
    if timeout <= 0:
        raise ValueError("NIST timeout must be positive.")
    if temp_root is None:
        temp_root = Path(tempfile.gettempdir())
        if not str(temp_root).isascii():
            temp_root = Path("C:/Temp")
    if not str(temp_root).isascii() or not temp_root.is_dir():
        raise ValueError("Supply an existing ASCII-only directory using --nist-temp.")
    msp.read_text(encoding="ascii")
    output.mkdir()
    lock = directory / "nist_cli.lock"
    previous_lock = directory / "gcms_ri_search.lock"
    if previous_lock.exists():
        raise ValueError(f"An earlier NIST bridge lock exists at {previous_lock}. Check NIST is idle before clearing it.")
    try:
        lock.touch(exist_ok=False)
    except FileExistsError as exc:
        raise ValueError(f"NIST bridge is locked at {lock}. Check NIST is idle before clearing it.") from exc
    stage = None
    started = False
    success = False
    originals = {}
    ready, result = directory / "SRCREADY.TXT", directory / "SRCRESLT.TXT"
    autoimp = directory / "AUTOIMP.MSD"
    try:
        stage = Path(tempfile.mkdtemp(prefix="gcms_ri_", dir=temp_root))
        for file in (autoimp, ready, result):
            originals[file] = file.read_bytes() if file.exists() else None
        backup = output / "previous_bridge_files"
        backup.mkdir()
        for file, data in originals.items():
            if data is not None:
                (backup / file.name).write_bytes(data)
        shutil.copyfile(directory / "nistms.INI", output / "settings_before.ini")
        shutil.copyfile(msp, stage / "input.msp")
        locator = stage / "locator.txt"
        locator.write_text(f"{stage / 'input.msp'} APPEND\n", encoding="ascii")
        ready.unlink(missing_ok=True)
        result.unlink(missing_ok=True)
        autoimp.write_text(f"{locator}\n", encoding="ascii")
        now = time.monotonic()
        proc = subprocess.Popen([str(executable), "/INSTRUMENT", "/PAR=2"], cwd=directory,
                                creationflags=subprocess.CREATE_NO_WINDOW)
        started = True
        while not ready.is_file():
            if time.monotonic() - now >= timeout:
                raise ValueError(f"NIST did not finish in {timeout:g} s. Its lock and staged input were retained.")
            time.sleep(0.25)
        data = result.read_bytes()
        if not data.strip():
            raise ValueError("NIST returned an empty report.")
        path = output / "results.txt"
        path.write_bytes(data)
        (output / "search.json").write_text(json.dumps({
            "program": str(executable), "flags": ["/INSTRUMENT", "/PAR=2"],
            "elapsed_seconds": round(time.monotonic() - now, 2),
            "launcher_exit_code": proc.poll(),
        }, indent=2), encoding="utf-8")
        success = True
        return path
    finally:
        try:
            for file, before in originals.items():
                if before is None:
                    file.unlink(missing_ok=True)
                else:
                    file.write_bytes(before)
        finally:
            if (success or not started) and stage is not None:
                shutil.rmtree(stage, ignore_errors=True)
            if success or not started:
                lock.unlink(missing_ok=True)
