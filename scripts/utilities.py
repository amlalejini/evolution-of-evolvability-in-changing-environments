import errno, os, csv

def mkdir_p(path):
    """
    This is functionally equivalent to the mkdir -p [fname] bash command
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def extract_params_cmd_log(path):
    """
    Extract Avida parameters from log of command used to run Avida.
    The log should contain only the text used to run Avida.

    e.g. a cmd log containing, './avida -set RANDOM_SEED 100' would return {"RANDOM_SEED": 100}
    """
    content = None
    with open(path, "r") as fp:
        content = fp.read().strip()
    content = content.replace("./avida", "")
    params = [param.strip() for param in content.split("-set") if param.strip() != ""]
    cfg = {param.split(" ")[0]:param.split(" ")[1] for param in params}
    return cfg

def read_avida_dat_file(path, backfill_missing_fields=False):
    """
    Parse an Avida .dat file (works with .spop files, too).
    Return parsed data as a list of dictionaries (where each dictionary is indexed by column name).
    """
    content = None
    with open(path, "r") as fp:
        content = fp.read().strip().split("\n")
    legend_start = 0
    legend_end = 0
    # Where does the legend table start?
    for line_i in range(0, len(content)):
        line = content[line_i].strip()
        if line == "# Legend:":         # Handles analyze mode detail files.
            legend_start = line_i + 1
            break
        if "#  1:" in line:             # Handles time.dat, mutation_rate.dat, etc.
            legend_start = line_i
            break
    # For each line in legend table, extract field
    fields = []
    for line_i in range(legend_start, len(content)):
        line = content[line_i].strip()
        if line == "":
            legend_end = line_i
            break
        # patch 3-input logic tasks because avida file format is nonsense
        if "Logic 3" in line:
            line = line.split("(")[0]

        fields.append( line.split(":")[-1].strip().lower().replace(" ", "_") )
    data = []
    for line_i in range(legend_end, len(content)):
        line = content[line_i].strip()
        if line == "": continue
        data_line = line.split(" ")
        if len(data_line) > len(fields):
            print("found more items than there are fields!")
            print(fields)
            print(data_line)
            exit(-1)
        elif backfill_missing_fields:
            num_backfill = len(fields) - len(data_line)
            for _ in range(num_backfill): data_line.append("")
        elif len(data_line) != len(fields):
            print("data fields mismatch!")
            print(fields)
            print(data_line)
            exit(-1)
        data.append({field:value for field,value in zip(fields, data_line)})
    return data

def read_avida_inst_set_file(inst_set_fpath):
    """
    Reads an avida instruction set file.
    Returns the mapping from instruction name to abbreviation (used in Avida genome sequence output).
    """
    content = None
    with open(inst_set_fpath, "r") as fp:
        content = fp.read().strip().split("\n")
    instructions = []
    char_lookup = []
    for line in content:
        line = line.strip()
        if line[:4] == "INST" and (not line.startswith("INSTSET")):
            instructions.append(line.split(" ")[1])
    for i in range(0, len(instructions)):
        if i < 26:
            char_lookup.append(chr(ord('a')+i))
        else:
            char_lookup.append(chr(ord('A')+(i-26)))
    return {inst:c for inst,c in zip(instructions, char_lookup)}

def read_csv(file_path):
    """
    Reads a csv, return contents as a list of dictionaries where each dictionary is a row in the csv that is indexed by column names.
    """
    content = None
    with open(file_path, "r") as fp:
        content = fp.read().strip().split("\n")
    header = content[0].split(",")
    content = content[1:]
    lines = [{header[i]: l[i] for i in range(len(header))} for l in csv.reader(content, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)]
    return lines

def simple_match_coeff(a, b):
    """
    Compute the simple matching coefficient between a and b.
    """
    if len(a) != len(b):
        print(f"Length mismatch! {a} {b}")
        exit(-1)
    return sum(ai==bi for ai,bi in zip(a,b))

def hamming_dist(a, b):
    """
    Compute the hamming distance between a and b.
    """
    if len(a) != len(b):
        print(f"Length mismatch! {a} {b}")
        exit(-1)
    return sum(ai!=bi for ai,bi in zip(a,b))