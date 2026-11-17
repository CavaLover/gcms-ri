using System;
using System.Globalization;
using System.IO;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Text;
using Agilent.MassSpectrometry.MSDChem;

internal static class AgilentReader
{
    private static string assemblyDirectory;

    private static Assembly ResolveAssembly(object sender, ResolveEventArgs args)
    {
        string path = Path.Combine(assemblyDirectory, new AssemblyName(args.Name).Name + ".dll");
        return File.Exists(path) ? Assembly.LoadFrom(path) : null;
    }

    private static int Main(string[] args)
    {
        if (args.Length != 3)
        {
            Console.Error.WriteLine("Usage: AgilentReader <assemblies> <data.ms> <scans.jsonl>");
            return 2;
        }
        try
        {
            assemblyDirectory = Path.GetFullPath(args[0]);
            AppDomain.CurrentDomain.AssemblyResolve += ResolveAssembly;
            Assembly.LoadFrom(Path.Combine(assemblyDirectory, "MSDChemDataFile.dll"));
            ReadScans(args[1], args[2]);
            return 0;
        }
        catch (Exception error)
        {
            Console.Error.WriteLine(error);
            return 1;
        }
    }

    private static void WriteArray(TextWriter writer, double[] values)
    {
        writer.Write('[');
        for (int i = 0; i < values.Length; i++)
        {
            if (i > 0) writer.Write(',');
            writer.Write(values[i].ToString("R", CultureInfo.InvariantCulture));
        }
        writer.Write(']');
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private static void ReadScans(string input, string output)
    {
        var file = new MSDChemDataMSFile();
        if (!file.OpenRead(input)) throw new IOException("Agilent could not open the input data.ms.");
        try
        {
            using (var writer = new StreamWriter(output, false, new UTF8Encoding(false)))
            {
                for (int i = 1; i <= file.NumberOfSpectra; i++)
                {
                    double rt;
                    ushort status;
                    double[] mz, intensity;
                    file.GetSpectrumByIndex1Based(i, out rt, out status, out mz, out intensity);
                    if (mz == null || intensity == null || mz.Length != intensity.Length)
                        throw new IOException("Invalid scan " + i.ToString(CultureInfo.InvariantCulture));
                    writer.Write("{\"scan\":" + i.ToString(CultureInfo.InvariantCulture)
                        + ",\"rt\":" + rt.ToString("R", CultureInfo.InvariantCulture) + ",\"mz\":");
                    WriteArray(writer, mz);
                    writer.Write(",\"intensity\":");
                    WriteArray(writer, intensity);
                    writer.WriteLine('}');
                }
            }
        }
        finally
        {
            file.Close();
        }
    }
}
