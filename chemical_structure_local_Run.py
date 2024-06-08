import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import io
from PIL import Image

# Function to retrieve SMILES string and compound details from PubChem
def get_compound_details(input_value):
    try:
        # Attempt to interpret input as a chemical compound name
        compounds = pcp.get_compounds(input_value, 'name')
        
        if compounds:
            # Extract SMILES string and compound details from the first result
            smiles = compounds[0].canonical_smiles
            compound = compounds[0]
            name = compound.iupac_name.replace(';', ' ') if compound.iupac_name else input_value
            return smiles, compound, name
        else:
            # If not found as a compound name, assume input is SMILES string
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(input_value))
            return smiles, None, None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None, None, None

# Function to visualize chemical structure based on input type
def visualize_chemical(input_value):
    smiles, compound, name = get_compound_details(input_value)
    
    if smiles:
        # Generate RDKit molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is not None:
            # Convert RDKit molecule to a PNG image
            img = Draw.MolToImage(mol, size=(500, 500))
            img_byte_arr = io.BytesIO()
            img.save(img_byte_arr, format='PNG')
            img_str = img_byte_arr.getvalue()

            # Display chemical information based on input type
            print("Chemical Information")
            print(f"- Name: {name}")
            if compound:
                print(f"- Formula: {compound.molecular_formula}")
                print(f"- Molecular Weight: {compound.molecular_weight}")

            # Display chemical structure image using Matplotlib
            img = Image.open(io.BytesIO(img_str))
            plt.imshow(img)
            plt.axis('off')
            plt.show()
        else:
            print("Failed to generate molecule from input.")
    else:
        print("Failed to retrieve SMILES for the input.")

# Main function to run the script
def main():
    print("Chemical Structure Plotter")
    print("Enter the name of the chemical compound or its SMILES string.")

    # Get user input for compound name or SMILES string
    input_value = input("Enter compound name or SMILES: ")

    # Visualize compound
    visualize_chemical(input_value)

if __name__ == "__main__":
    main()
