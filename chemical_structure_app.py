import streamlit as st
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Draw
import io

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
        st.error(f"An error occurred: {e}")
        return None, None, None

# Function to visualize chemical structure based on input type
def visualize_chemical(input_value):
    smiles, compound, name = get_compound_details(input_value)
    
    if smiles:
        # Generate RDKit molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is not None:
            # Convert RDKit molecule to a PNG image
            img = Draw.MolToImage(mol, size=(800, 500))
            img_byte_arr = io.BytesIO()
            img.save(img_byte_arr, format='PNG')
            img_str = img_byte_arr.getvalue()

            # Display chemical information based on input type
            st.subheader("Chemical Information")
            st.markdown(f"- **Name:** {name}")
            st.markdown(f"- **Formula:** {compound.molecular_formula}")
            st.markdown(f"- **Molecular Weight:** {compound.molecular_weight}")

            # Display chemical structure image
            st.subheader("Chemical Structure")
            st.image(img_str, caption='Chemical Structure', use_column_width=True)
        else:
            st.error("Failed to generate molecule from input.")
    else:
        st.error("Failed to retrieve SMILES for the input.")

# Main function to create Streamlit app
def main():
    # Set page title and favicon
    st.set_page_config(page_title="Chemical Structure Plotter", page_icon=":microscope:")

    # Set app layout and appearance
    st.title("Chemical Structure Plotter")
    st.sidebar.image("https://upload.wikimedia.org/wikipedia/commons/thumb/f/f0/Nucleus_drawing.svg/800px-Nucleus_drawing.svg.png", width=100)
    st.sidebar.markdown("**Welcome to the Chemical Structure Plotter!**")
    st.sidebar.markdown("Enter the name of the chemical compound or its SMILES string in the text box below.")
    st.sidebar.markdown("Click the 'Visualize' button to see its chemical structure and details.")
    st.sidebar.markdown("---")
    st.sidebar.markdown("Follow me on:")
    st.sidebar.markdown(
        '<a href="https://github.com/ariktheone"><img src="https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcSmrGmeBv3SOLSKz6OlTVlVYkfH9_W3BBgdrA&s" width="40" height="40"></a>'
        '<a href="https://www.linkedin.com/in/arijitmondal30/"><img src="https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcRokEYt0yyh6uNDKL8uksVLlhZ35laKNQgZ9g&s" width="40" height="40"></a>'
        '<a href="https://twitter.com/ArijitM63768876"><img src="https://cdn-icons-png.flaticon.com/512/124/124021.png" width="40" height="40"></a>',
        unsafe_allow_html=True
    )

    # Get user input for compound name or SMILES string
    input_value = st.text_input("Enter compound name or SMILES")

    # Visualize compound on button click
    if st.button("Visualize"):
        visualize_chemical(input_value)

if __name__ == "__main__":
    main()

