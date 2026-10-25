import streamlit as st 
import pandas as pd 
from PIL import Image 
import subprocess 
import os 
import pickle
from pathlib import Path 
import base64
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

BASE_DIR = Path(__file__).resolve().parent

INPUT_FILE = BASE_DIR / "molecule.smi"
OUTPUT_FILE = BASE_DIR / "descriptor_output.csv"

PADEL_JAR = Path(
    r"D:\ANDINA\Course\Bioinformatics\Python\EGFR Drug Discovery\downloaded_repo\bioactivity-prediction-app-main\PaDEL-Descriptor\PaDEL-Descriptor.jar"
)

FINGERPRINT_XML = Path(
    r"D:\ANDINA\Course\Bioinformatics\Python\EGFR Drug Discovery\fingerprints_xml\PubchemFingerprinter.xml"
)

# Molecular descriptor calculator
def desc_calc():
    command = [
        "java",
        "-Xms2G",
        "-Xmx2G",
        "-Djava.awt.headless=true",
        "-jar", str(PADEL_JAR),
        "-removesalt",
        "-standardizenitro",
        "-fingerprints",
        "-descriptortypes", str(FINGERPRINT_XML),
        "-dir", str(BASE_DIR),
        "-file", str(OUTPUT_FILE),
    ]

    # Check that the required PaDEL files exist
    if not PADEL_JAR.exists():
        raise FileNotFoundError(f"PaDEL JAR not found: {PADEL_JAR}")

    if not FINGERPRINT_XML.exists():
        raise FileNotFoundError(
            f"Fingerprint XML not found: {FINGERPRINT_XML}"
        )

    # Run PaDEL and capture any error messages
    result = subprocess.run(
        command,
        cwd=BASE_DIR,
        capture_output=True,
        text=True
    )

    # Stop if PaDEL failed
    if result.returncode != 0:
        st.error("PaDEL descriptor calculation failed.")
        st.code(result.stderr or result.stdout)
        st.stop()

    # Confirm that the output file was created
    if not OUTPUT_FILE.exists():
        st.error("PaDEL finished but did not create descriptors_output.csv.")
        st.code(result.stdout)
        st.stop()

    # Remove the temporary input only after successful calculation
    if INPUT_FILE.exists():
        INPUT_FILE.unlink()


# Lipinski characteristics calculator
def calculate_lipinski(smiles):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Invalid SMILES string.")

    return {
        "Molecular Weight": Descriptors.MolWt(mol),
        "LogP": Crippen.MolLogP(mol),
        "H-bond Donors": Lipinski.NumHDonors(mol),
        "H-bond Acceptors": Lipinski.NumHAcceptors(mol),
    }

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building
def build_model(input_data, load_data):
    # Load the trained EGFR regression model
    with open(BASE_DIR / "EGFR_model.pkl", "rb") as f:
        model = pickle.load(f)

    # Predict pIC50 using the PaDEL features
    predictions = model.predict(input_data)

    # Calculate Lipinski descriptors for each input SMILES
    lipinski_results = []

    for smiles in load_data["SMILES"]:
        properties = calculate_lipinski(smiles)
        lipinski_results.append(properties)

    lipinski_df = pd.DataFrame(lipinski_results)

    # Combine molecule IDs, predictions, and descriptors
    results = pd.DataFrame({
        "ChEMBL ID": load_data["ChEMBL ID"].values,
        "Predicted pIC50": predictions
    })

    results = pd.concat(
        [results, lipinski_df],
        axis=1
    )
    return results

# Logo image
image = Image.open(
    r"D:\ANDINA\Course\Bioinformatics\Python\EGFR Drug Discovery\Drug Discovery App\logo.png"
)

st.image(image, width="stretch")

# Page title
st.markdown("""
# Bioactivity Prediction App (EGFR - Epidermal Growth Factor Receptor)

This app allows you to predict the bioactivity towards inhibting the `EGFR` enzyme.
EGFR is a receptor tyrosine kinase involved in cell growth and signalling. 
Its dysregulation or mutation is relevant to several cancers and many small-molecule inhibitors have been studied against it.

**Credits**
- App built in `Python` + `Streamlit` by [Chanin Nantasenamat](https://medium.com/@chanin.nantasenamat) (aka [Data Professor](http://youtube.com/dataprofessor))
- Descriptor calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Read the Paper]](https://doi.org/10.1002/jcc.21707).
---
""")

# Sidebar
with st.sidebar.header('1. Upload your CSV data'):
    uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])
    st.sidebar.markdown("""
[Example input file](https://raw.githubusercontent.com/dataprofessor/bioactivity-prediction-app/main/example_acetylcholinesterase.txt)
""")

if st.sidebar.button("Predict"):
    if uploaded_file is None:
        st.warning("Please upload a TXT file first.")
        st.stop()

    # Read the uploaded SMILES and ChEMBL IDs
    load_data = pd.read_csv(
        uploaded_file,
        sep=r"\s+",
        header=None,
        names=["SMILES", "ChEMBL ID"],
        dtype=str
    )

    if load_data.empty:
        st.error("The uploaded file is empty.")
        st.stop()

    # Display original input
    st.header("Original Input Data")
    st.dataframe(load_data)

    # Write input file for PaDEL
    load_data[["SMILES", "ChEMBL ID"]].to_csv(
        INPUT_FILE,
        sep="\t",
        header=False,
        index=False
    )

    # Calculate PaDEL fingerprints
    with st.spinner("Calculating molecular fingerprints..."):
        desc_calc()

    if not OUTPUT_FILE.exists():
        st.error("PaDEL output file was not created.")
        st.stop()

    desc = pd.read_csv(OUTPUT_FILE)

    st.header("Calculated Molecular Descriptors")
    st.dataframe(desc)
    st.write("Descriptor matrix shape:", desc.shape)

    # Select the exact features used during model training
    Xlist = pd.read_csv(
        BASE_DIR / "descriptor_list.csv"
    ).columns.tolist()

    missing_columns = [
        column for column in Xlist
        if column not in desc.columns
    ]

    if missing_columns:
        st.error(
            f"Missing model features: {missing_columns[:10]}"
        )
        st.stop()

    desc_subset = desc[Xlist]

    # Make predictions and calculate Lipinski properties
    results = build_model(desc_subset, load_data)

    st.header("EGFR Bioactivity Prediction")
    st.dataframe(results, use_container_width=True)

    # Download the complete result table
    st.download_button(
        label="Download Predictions and Lipinski Descriptors",
        data=results.to_csv(index=False).encode("utf-8"),
        file_name="EGFR_predictions_with_lipinski.csv",
        mime="text/csv"
    )

else:
    st.info("Upload a TXT file containing SMILES and ChEMBL IDs.")

# to run it, paste this on terminal line
# & "C:\Users\Andina Angelia\AppData\Local\Python\pythoncore-3.14-64\python.exe" -m streamlit run "D:\ANDINA\Course\Bioinformatics\Python\EGFR Drug Discovery\Drug Discovery App\App.py"