import os
import re
import shutil

# ----------------- PARSE GAUSSIAN OUT -----------------

def parse_gaussian_out(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()
    # --- CLEAN TEXT FOR BROKEN GAUSSIAN LINES ---
    # Remove line breaks inside parentheses
    # --- CHECK NORMAL TERMINATION ---
    terminated = bool(re.search(r"Normal termination of Gaussian", text))
    if not terminated:
        print(f"WARNING: Gaussian did not terminate properly: {path}")
        exit()
    text = re.sub(r"\(\s*\n\s*", "(", text)
    text = re.sub(r"\s*\n\s*\)", ")", text)

    # Remove mid‑word line breaks (e.g. TetraHy\n droFuran)
    text = re.sub(r"([A-Za-z])\s*\n\s*([A-Za-z])", r"\1\2", text)

    # Collapse multiple spaces
    text = re.sub(r"\s+", " ", text)

    # --- FUNCTIONAL ---
    m = re.search(r"(?i)#p?\s+([A-Za-z0-9-]+)", text)
    functional = m.group(1).upper() if m else "UNKNOWN"

    if functional == "PBE1PBE":
        functional = "PBE0"
    if functional == "CAM-B3LYP":
        functional = "CAMB3LYP"
    if functional in ["WB97XD", "WB97X-D"]:
        functional = "wB97XD"

    # --- SOLVENT ---
    if re.search(r"(?i)SCRF=\(Solvent\s*=\s*TetraHydroFuran\)", text):
        solvent = "THF"
    elif re.search(r"(?i)SCRF=\(Solvent\s*=\s*Dichloromethane\)", text):
        solvent = "DCM"
    else:
        solvent = "Vacuum"

    # --- STATE ---

    # TDA (singlets or triplets)
    if re.search(r"(?i)\btda\s*(=)?\s*\(", text):
        state = "TDA"

    # TD singlets → S1
    elif re.search(r"(?i)\btd\s*(=)?\s*\(\s*singlets", text):
        state = "S1"

    # TD (triplets or unspecified)
    elif re.search(r"(?i)\btd\s*(=)?\s*\(", text):
        state = "TD"

    # Ground-state singlet
    elif re.search(r"(?i)Multiplicity\s*=\s*1", text):
        state = "S0"

    # Ground-state triplet
    elif re.search(r"(?i)Multiplicity\s*=\s*3", text):
        state = "TDFT"

    else:
        state = "UNKNOWN"

    return [solvent, state, functional]


# ----------------- NORMALIZATION / ALIASES -----------------

FUNCTIONALS = ["B3LYP", "PBE0", "TPSSH", "M062X", "CAMB3LYP", "wB97XD"]
SOLVENTS    = ["Vacuum", "DCM", "THF"]
STATES      = ["S0", "TDFT", "TDA", "TD", "S1"]
ISOMERS     = ["Fac", "Mer", "Cis", "Trans"]

functional_to_index = {f: i for i, f in enumerate(FUNCTIONALS)}
solvent_to_index    = {s: i for i, s in enumerate(SOLVENTS)}
state_to_index      = {s: i for i, s in enumerate(STATES)}
isomer_to_index     = {i: j for j, i in enumerate(ISOMERS)}

ALIASES = {
    # Functionals
    "PBE1": "PBE0",
    "PBE": "PBE0",
    "CAM-B3LYP": "CAMB3LYP",
    "WB97XD": "wB97XD",
    "WB97X-D": "wB97XD",
    "TPPSh": "TPSSH",
    "TPSSh" : "TPSSH",
    "TPPsh" : "TPSSH",

    # Solvents
    "GP": "Vacuum",
    "GasPhase": "Vacuum",
    "SCM-DCM": "DCM",
    "PCM-DCM": "DCM",
    "SCM-THF": "THF",
    "PCM-THF": "THF",

    # States
    "T-DFT": "TDFT",
    "T-TDA": "TDA",
    "S1toS0": "S1",
    "T1TDA": "TDA",
    "T1TD": "TD",
    "T": "TDFT",
    "T1-TDA" : "TDA",
    "T1-TD" : "TD",
    "T-TD" : "TD",

    # Isomers
    "fac": "Fac",
    "mer": "Mer",
    "cis": "Cis",
    "trans": "Trans",
}

def normalize_token(token):
    return ALIASES.get(token, token)


# ----------------- FILENAME PARSING -----------------

def parse_filename(fname):
    base = os.path.splitext(os.path.basename(fname))[0]
    parts = base.split("_")

    solvent = None
    state = None
    isomer = None
    functional = None

    for p in parts:
        p_norm = normalize_token(p)

        if p_norm in FUNCTIONALS:
            functional = p_norm
        elif p_norm in SOLVENTS:
            solvent = p_norm
        elif p_norm in STATES:
            state = p_norm
        elif p_norm in ISOMERS:
            isomer = p_norm

    return {
        "filename": fname,
        "solvent": solvent,
        "state": state,
        "isomer": isomer,
        "functional": functional
    }


# ----------------- CONSISTENCY CHECK -----------------

def check_consistency(filename_info, out_info):
    """
    filename_info: dict from parse_filename()
    out_info: [solvent, state, functional] from parse_gaussian_out()
    """

    out_solvent, out_state, out_functional = out_info
    warnings = []

    # --- Solvent ---
    if filename_info["solvent"] is None:
        warnings.append("Filename missing solvent information")
    elif filename_info["solvent"] != out_solvent:
        warnings.append(f"Solvent mismatch: filename={filename_info['solvent']} vs out={out_solvent}")

    # --- State ---
    if filename_info["state"] is None:
        warnings.append("Filename missing state information")
    elif filename_info["state"] != out_state:
        warnings.append(f"State mismatch: filename={filename_info['state']} vs out={out_state}")

    # --- Functional ---
    if filename_info["functional"] is None:
        warnings.append("Filename missing functional information")
    elif filename_info["functional"] != out_functional:
        warnings.append(f"Functional mismatch: filename={filename_info['functional']} vs out={out_functional}")

    return warnings



# ----------------- MAIN -----------------

if __name__ == "__main__":
    input_folder = os.path.join(os.getcwd(), "Input", "out_to_info", "raw", "Mer")
    input_folder = r"C:\Users\leosa\Desktop\Internship\All out\Old_name\Mer"
    files = [
        os.path.join(input_folder, f)
        for f in os.listdir(input_folder)
        if f.endswith(".out")
    ]
    warning = 0
    for f in files:
        filename_info = parse_filename(os.path.basename(f))
        out_info = parse_gaussian_out(f)

        print("\nFILE:", os.path.basename(f))
        print("Filename info:", filename_info)
        print("Out info     :", out_info)
        warnings = check_consistency(filename_info, out_info)

        if warnings:
            print("⚠️  CONSISTENCY WARNINGS:")
            warning += 1
            for w in warnings:
                print("   -", w)
        else:
            print("✓ Filename and .out are consistent")
            shutil.copy(f, os.path.join(os.getcwd(), "Input", "out_to_info", "Renamed", filename_info["isomer"],f"{filename_info["isomer"]}_{out_info[0]}_{out_info[1]}_{out_info[2]}.out"))
print(f"Total warnings: {warning}")