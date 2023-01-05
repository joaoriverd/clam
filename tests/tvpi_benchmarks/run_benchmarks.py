import subprocess
import os
import csv

# Paths to benchmarks and name
benchmarks = [
              # ["name", "rel_path/benchmark_file.c"],
              ["artificial", "papers/artificial_benchmarks.c"],
              ["fppoly", "papers/fppoly_benchmarks.c"],
              ["fptylor_paper", "papers/fptylor_benchmarks.c"],
              ["raicp", "papers/raicp_benchmarks.c"],
              ["rosa", "fpbench/rosa_benchmarks.c"],
              ["daisy", "fpbench/daisy_benchmarks.c"],
              ["fptaylor", "fpbench/fptaylor_benchmarks.c"],
             ]

clam_exec = "/home/joao/Documents/repo/static_analysis/clam/install/bin/clam.py"
benchmarks_path = "/home/joao/Documents/repo/static_analysis/clam/tests/tvpi_benchmarks/"
results_path = "/home/joao/Documents/repo/static_analysis/clam/tests/tvpi_benchmarks/results/"
analysis_path = "/home/joao/Documents/repo/static_analysis/clam/tests/tvpi_benchmarks/results/analysis/"
analysis_out_prefix = analysis_path + "/analysis_out_"
analysis_err_prefix = analysis_path + "/analysis_err_"
clam_attr = ["--crab-verbose=5",
             "--crab-widening-delay=1",
             "--crab-narrowing-iterations=2",
             "--crab-dom=tvpi",
             "--crab-check=assert",
             "--save-temps",
             "--temp-dir=" + analysis_path,
             # "--llvm-dot-cfg",
             ]

header_csv = ["function", "var", "lower_bound", "upper_bound", "range"]


def run_command(cmd):
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return out.decode("utf-8"), err.decode("utf-8")

def substring_after(s, delim):
    return s.partition(delim)[2]

def substring_before(s, delim):
    return s.partition(delim)[0]

def substring_between(s, delim1, delim2):
    t = substring_after(s, delim1)
    return substring_before(t, delim2)

def write_to_file(s, file):
    text_file = open(file, "w")
    n = text_file.write(s)
    text_file.close()

if __name__ == "__main__":
    for b in benchmarks:
        # Prepare command to run benchmark
        benchmark_name = b[0]
        benchmark_file = benchmarks_path + b[1]
        cmd = [clam_exec, benchmark_file] + clam_attr + ["-ocrab=" + analysis_path + benchmark_name]
        clam_out, clam_err = run_command(cmd)
        write_to_file(clam_out, analysis_out_prefix + benchmark_name)
        write_to_file(clam_err, analysis_err_prefix + benchmark_name)
        print("\nRunning Benchmark: " + benchmark_name)

        # Search for variables where crab_range is calculated within the analysis output and
        # save such variables in a csv file.
        f = open(results_path + benchmark_name + ".csv", "w")
        if not f:
            print("An error occurred...")
        writer = csv.writer(f)
        writer.writerow(header_csv)

        clam_out = clam_out.splitlines()
        current_function = ""
        for l in clam_out:
            # Check if the function name is found
            function_name = substring_after(l, "Starting CFG construction for ")
            if function_name:
                print(function_name)
                current_function = function_name
                continue

            crab_range = substring_after(l, "crab_range: ")
            if crab_range:
                print(" " + crab_range)
                # Prepare data for csv row
                var_name = substring_before(crab_range, " : ")
                range_lower = substring_between(crab_range, "[", ",")
                range_upper = substring_between(crab_range, ", ", "]")
                range_complete = substring_after(crab_range, " : ")
                csv_row = [current_function, var_name, range_lower, range_upper, range_complete]
                writer.writerow(csv_row)
