"""Use a user-supplied Agilent installation to read .D/data.ms on Windows."""

import os
from importlib.resources import as_file, files
from pathlib import Path
import shutil
import subprocess
import tempfile


def compile_reader(directory, temp, arch):
    if os.name != "nt":
        raise ValueError("Native Agilent reading requires Windows. Import an exported peak table and MSP instead.")
    library = directory / "MSDChemDataFile.dll"
    if not library.is_file():
        raise ValueError(f"MSDChemDataFile.dll was not found in {directory}.")
    framework = "Framework64" if arch == "x64" else "Framework"
    compiler = Path(os.environ.get("SystemRoot", "C:/Windows")) / "Microsoft.NET" / framework / "v4.0.30319" / "csc.exe"
    if not compiler.is_file():
        raise ValueError(f"The .NET Framework C# compiler was not found at {compiler}.")
    reader = temp / "AgilentReader.exe"
    with as_file(files("gcms_ri").joinpath("AgilentReader.cs")) as source:
        _run([str(compiler), "/nologo", f"/platform:{arch}", f"/out:{reader}",
              f"/reference:{library}", str(source)], cwd=temp, timeout=60)
    return reader


def _run(command, cwd, timeout, env=None):
    result = subprocess.run(command, cwd=cwd, capture_output=True, timeout=timeout, env=env,
                            creationflags=subprocess.CREATE_NO_WINDOW)
    if result.returncode:
        message = b"\n".join([result.stdout, result.stderr]).decode("mbcs", errors="replace").strip()
        raise ValueError(f"Agilent reader exited with status {result.returncode}: {message}")


def native_extract(samples, directory, work, arch, timeout):
    """Read temporary copies so the vendor reader never opens the untouched source."""
    env = dict(os.environ)
    env["PATH"] = str(directory) + os.pathsep + env.get("PATH", "")
    with tempfile.TemporaryDirectory(prefix="gcms-ri-reader-") as temp_name:
        temp = Path(temp_name)
        reader = compile_reader(directory, temp, arch)
        for label, sample in samples:
            original = sample / "data.ms" if sample.is_dir() else sample
            if not original.is_file() or original.suffix.lower() != ".ms":
                raise ValueError(f"Missing Agilent .D/data.ms: {sample}")
            copy = temp / "data.ms"
            shutil.copyfile(original, copy)
            partial = work / f"{label}.scans.partial"
            final = work / f"{label}.scans.jsonl"
            _run([str(reader), str(directory), str(copy), str(partial)], cwd=temp, timeout=timeout, env=env)
            partial.replace(final)
            yield label, sample, final
