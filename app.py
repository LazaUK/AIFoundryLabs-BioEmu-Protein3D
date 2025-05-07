# BioEmu AI Model Demo Web App - Simplified Version
# ----------------------------------------------
import gradio as gr
import os
import tempfile
import shutil
import glob
from pathlib import Path
import logging
import html  # For escaping HTML for iframe srcdoc

# Attempt to import BioEmu and handle potential errors
try:
    from bioemu.sample import main as bioemu_sample
except ImportError as e:
    logging.exception("Failed to import bioemu.sample - Please ensure BioEmu is installed correctly.")
    print(f"Error importing bioemu: {e}")
    print("Please ensure the bioemu package and its dependencies are installed correctly in your environment.")
    raise  # Stop execution if bioemu is essential

# Import other necessary libraries
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import py3Dmol # For 3D visualization
import json
import numpy as np
import traceback

# Suppress Biopython PDB warnings for cleaner output
warnings.filterwarnings("ignore", category=PDBConstructionWarning)

# --- Main Application Class ---
class BioEmuDemo:
    """Handles BioEmu execution and temporary file management."""
    def __init__(self):
        """Initialize the BioEmu demo interface."""
        self.temp_dirs = []
        logging.info("BioEmuDemo initialized.")

    def __del__(self):
        """Clean up temporary directories on exit."""
        logging.info(f"Cleaning up {len(self.temp_dirs)} temporary directories...")
        for temp_dir in self.temp_dirs:
            if os.path.exists(temp_dir):
                try:
                    shutil.rmtree(temp_dir)
                    logging.info(f"Successfully removed temp directory: {temp_dir}")
                except Exception as e:
                    logging.error(f"Error cleaning up temp directory {temp_dir}: {e}")
            else:
                logging.warning(f"Attempted to clean up non-existent directory: {temp_dir}")

    def run_sample(self, sequence, num_samples=5):
        """
        Runs BioEmu sampling, processes results and generates output for Gradio.
        (This method is called by the run_and_update_wrapper which handles progress)
        """
        iframe_html = None
        score_plot = None
        output_dir = None

        if not sequence or len(sequence.strip()) == 0:
            logging.warning("run_sample called with empty sequence.")
            return "Please enter a protein sequence.", iframe_html, score_plot, output_dir

        original_sequence_head = sequence[:60].replace('\n', ' ')
        if sequence.startswith('>'):
            lines = sequence.strip().split('\n')
            sequence = ''.join(lines[1:]).replace(' ', '').replace('\t', '').replace('\r', '')
        else:
            sequence = sequence.strip().replace(' ', '').replace('\t', '').replace('\r', '')
        logging.info(f"Processing sequence (len={len(sequence)}): {sequence[:60]}...")

        valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
        upper_sequence = sequence.upper()
        if not all(aa in valid_aa for aa in upper_sequence):
            invalid_chars = sorted(list(set(aa for aa in upper_sequence if aa not in valid_aa)))
            err_msg = f"Invalid amino acids found: {', '.join(invalid_chars)}"
            logging.error(f"{err_msg} in sequence starting with: {original_sequence_head}")
            return err_msg, iframe_html, score_plot, output_dir

        try:
            output_dir = tempfile.mkdtemp(prefix="bioemu_")
            self.temp_dirs.append(output_dir)
            logging.info(f"Created temp output directory: {output_dir}")

            logging.info(f"Running BioEmu for sequence: {sequence[:20]}..., N={num_samples}")
            bioemu_sample(
                sequence=sequence,
                num_samples=int(num_samples), # Ensure num_samples is int
                output_dir=output_dir,
                filter_samples=False
            )
            logging.info("BioEmu filter_samples set to False.")
            logging.info(f"BioEmu process finished. Checking output in {output_dir}")

            display_pdb_path = None
            num_generated = 0
            topology_pdb = os.path.join(output_dir, "topology.pdb")
            individual_pdbs = sorted(glob.glob(os.path.join(output_dir, "rank*.pdb")))

            if individual_pdbs:
                 try:
                     individual_pdbs.sort(key=lambda f: int(os.path.basename(f).split('.')[0].replace('rank','')))
                 except ValueError:
                      logging.warning(f"Could not sort ranked PDBs by number, using default glob order: {[os.path.basename(f) for f in individual_pdbs]}")
                 display_pdb_path = individual_pdbs[0]
                 num_generated = len(individual_pdbs)
                 logging.info(f"Found {num_generated} ranked PDB(s). Using: {os.path.basename(display_pdb_path)}")
            elif os.path.exists(topology_pdb):
                 display_pdb_path = topology_pdb
                 logging.info(f"No ranked PDBs found. Using topology.pdb.")
                 npz_files = glob.glob(os.path.join(output_dir, "batch_*.npz"))
                 if npz_files:
                     try: num_generated = sum(np.load(f)['pos'].shape[0] for f in npz_files)
                     except Exception as e:
                         logging.warning(f"Could not accurately count samples from npz files: {e}")
                         num_generated = int(num_samples)
                 else: num_generated = 1 # If topology.pdb exists, assume at least 1 sample from it
            else:
                 npz_files_exist = bool(glob.glob(os.path.join(output_dir, "*.npz")))
                 msg = "No displayable PDB file (topology.pdb or rank*.pdb) was found."
                 if npz_files_exist: msg += " Intermediate (.npz) files exist."
                 else: msg += " No .npz files found either."
                 logging.error(f"{msg} Path: {output_dir}")
                 return msg, iframe_html, score_plot, output_dir

            pdb_data = ""
            try:
                with open(display_pdb_path, 'r') as f: pdb_data = f.read()
                logging.info(f"Read PDB data (size={len(pdb_data)}) from {display_pdb_path}")
                if not pdb_data: raise ValueError("PDB file is empty.")
            except Exception as e:
                 err_msg = f"Error reading PDB file {display_pdb_path}: {e}"
                 logging.exception(err_msg)
                 return err_msg, iframe_html, score_plot, output_dir

            try:
                logging.info("Creating py3Dmol view...")
                viewer_html_content = self._create_molecule_viewer(pdb_data)

                if viewer_html_content and not viewer_html_content.startswith("<p>Failed"):
                    logging.info("Successfully generated py3Dmol HTML content.")
                    escaped_html = html.escape(viewer_html_content)
                    iframe_html = f'<iframe srcdoc="{escaped_html}" width="610px" height="510px" style="border: none;" title="3D Protein Structure Viewer"></iframe>'
                    logging.info("Wrapped py3Dmol HTML in iframe srcdoc.")
                elif viewer_html_content:
                    logging.error(f"py3Dmol view creation failed internally, returned error HTML.")
                    iframe_html = viewer_html_content
                else:
                    logging.error("py3Dmol _create_molecule_viewer returned None or empty string.")
                    iframe_html = "<p>Error: Failed to generate 3D viewer HTML content.</p>"

            except Exception as e:
                 err_msg = f"Error processing py3Dmol output or creating iframe: {e}"
                 logging.exception(err_msg)
                 iframe_html = f"<p>Error creating molecule viewer display: {e}</p>"

            result_summary = f"Generated {num_generated} structure(s) for sequence: {sequence}\n\n"
            result_summary += f"Structures saved in: {output_dir}\n"
            result_summary += f"Displaying structure: {os.path.basename(display_pdb_path)}"
            xtc_file = os.path.join(output_dir, "samples.xtc")
            if os.path.exists(topology_pdb) and os.path.exists(xtc_file):
                result_summary += f"\n\nFound topology.pdb and samples.xtc for trajectory analysis."
            else:
                 missing = []
                 if not os.path.exists(topology_pdb): missing.append("topology.pdb")
                 if not os.path.exists(xtc_file): missing.append("samples.xtc")
                 if missing: result_summary += f"\n\nNote: Did not find expected file(s): {', '.join(missing)}."

            logging.info("run_sample finished successfully.")
            return result_summary, iframe_html, score_plot, output_dir

        except Exception as e:
            logging.exception(f"Critical error during run_sample for sequence: {sequence[:60]}...")
            current_output_dir = locals().get('output_dir', None)
            return f"Unexpected Error during processing: {str(e)}", None, None, current_output_dir

    
    def _create_molecule_viewer(self, pdb_data):
        """Generates the raw HTML+JS string for the 3D viewer using py3Dmol."""
        try:
            logging.info("Attempting to create py3Dmol view...")
            viewer_width = 600
            viewer_height = 500
            view = py3Dmol.view(width=viewer_width, height=viewer_height)
            logging.info(f"py3Dmol.view object created, width={viewer_width}, height={viewer_height}.")

            view.addModel(pdb_data, "pdb")
            logging.info("PDB data added to model.")

            # Stick representation, colored by element
            view.setStyle({'stick': {'colorscheme': 'element'}})
            logging.info("Stick style applied.")

            # Set background and camera
            view.zoomTo()
            view.setBackgroundColor('white')
            
            # Add spin control
            view.spin(True)
            logging.info("Added spin control to view.")
            
            # Create a color legend for elements
            legend_html = """
            <div style="position: absolute; top: 10px; right: 10px; background-color: rgba(255,255,255,0.7); 
                        padding: 10px; border-radius: 5px; border: 1px solid #ccc; font-size: 12px;">
                <div style="font-weight: bold; margin-bottom: 5px;">Element Colors:</div>
                <div><span style="display: inline-block; width: 15px; height: 15px; background-color: #C8C8C8; vertical-align: middle;"></span> Carbon (C)</div>
                <div><span style="display: inline-block; width: 15px; height: 15px; background-color: #FF0000; vertical-align: middle;"></span> Oxygen (O)</div>
                <div><span style="display: inline-block; width: 15px; height: 15px; background-color: #8F8FFF; vertical-align: middle;"></span> Nitrogen (N)</div>
                <div><span style="display: inline-block; width: 15px; height: 15px; background-color: #FFC832; vertical-align: middle;"></span> Sulfur (S)</div>
                <div><span style="display: inline-block; width: 15px; height: 15px; background-color: #FFFFFF; vertical-align: middle; border: 1px solid #ccc;"></span> Hydrogen (H)</div>
            </div>
            """
            
            # Generate viewer HTML and add legend
            html_output = view._make_html()
            
            # Insert legend before the closing body tag
            if "</body>" in html_output:
                html_output = html_output.replace("</body>", f"{legend_html}</body>")
            else:
                html_output += legend_html
                
            logging.info("Generated HTML with element color legend.")
            return html_output
            
        except Exception as e:
            logging.exception("Error occurred inside _create_molecule_viewer")
            return f"<p>Failed to create 3D viewer internally: {html.escape(str(e))}</p>"

    def download_structures(self, output_dir):
        """Creates a zip archive of the generated structure files and logs."""
        if not output_dir or not os.path.exists(output_dir):
             gr.Warning("Output directory not found or invalid. Please generate structures first.")
             logging.warning(f"download_structures called with invalid/missing output_dir: {output_dir}")
             return None

        # Define zip path outside the staging directory to avoid zipping itself
        final_zip_dir = os.path.dirname(output_dir) if os.path.dirname(output_dir) and os.path.exists(os.path.dirname(output_dir)) else tempfile.gettempdir()
        final_zip_basename = os.path.join(final_zip_dir, f"bioemu_results_{os.path.basename(output_dir)}")
        final_zip_path = final_zip_basename + ".zip"

        logging.info(f"Attempting to create zip file: {final_zip_path} from dir: {output_dir}")

        files_to_zip = []
        patterns_to_check = ["rank*.pdb", "topology.pdb", "samples.xtc", "sequence.fasta", "*.npz", "*.log"]
        for pattern in patterns_to_check:
            found_files = glob.glob(os.path.join(output_dir, pattern))
            abs_found_files = [os.path.abspath(f) for f in found_files]
            abs_files_already_added = [os.path.abspath(f_existing) for f_existing in files_to_zip]
            new_files_to_add = [f for f_abs, f in zip(abs_found_files, found_files) if f_abs not in abs_files_already_added]
            if new_files_to_add:
                files_to_zip.extend(new_files_to_add)

        if not files_to_zip:
             gr.Warning("No relevant files found to zip in the output directory.")
             logging.warning(f"No files found to zip in {output_dir}.")
             return None

        staging_dir = None
        try:
            parent_temp_dir = tempfile.gettempdir()
            staging_dir = tempfile.mkdtemp(prefix="bioemu_zip_staging_", dir=parent_temp_dir)
            logging.info(f"Using staging directory for zip: {staging_dir}")

            copied_files_count = 0
            for f_path in files_to_zip:
                if os.path.exists(f_path):
                     try:
                         shutil.copy(f_path, os.path.join(staging_dir, os.path.basename(f_path)))
                         copied_files_count += 1
                     except Exception as copy_e: logging.warning(f"Could not copy file {os.path.basename(f_path)} to staging: {copy_e}")
                else: logging.warning(f"File {os.path.basename(f_path)} not found during copy for zipping.")

            if copied_files_count == 0:
                 gr.Warning("No files were successfully copied for zipping.")
                 logging.warning(f"No files copied to staging dir {staging_dir} from {output_dir}")
                 return None

            logging.info(f"Copied {copied_files_count} files to staging dir {staging_dir}.")

            if os.path.exists(final_zip_path):
                try:
                    os.remove(final_zip_path)
                    logging.info(f"Removed existing zip file: {final_zip_path}")
                except Exception as e_rm:
                    logging.warning(f"Could not remove existing zip file {final_zip_path}: {e_rm}. Attempting to overwrite.")


            shutil.make_archive(final_zip_basename, 'zip', root_dir=staging_dir) # Use root_dir for cleaner zip contents
            logging.info(f"Successfully created archive: {final_zip_path}")

            if os.path.exists(final_zip_path): return final_zip_path
            else:
                gr.Error("Failed to create zip file after archive operation (file not found).")
                logging.error(f"Zip file {final_zip_path} not found after make_archive.")
                return None
        except Exception as e:
            gr.Error(f"An error occurred while creating the zip file: {e}")
            logging.exception("Error during zipping process.")
            return None
        finally:
             if staging_dir and os.path.exists(staging_dir):
                 try:
                    shutil.rmtree(staging_dir)
                    logging.info(f"Successfully removed staging directory: {staging_dir}")
                 except Exception as cleanup_e: logging.error(f"Error cleaning staging dir {staging_dir}: {cleanup_e}")

# --- Sample Sequences Definition ---
SAMPLE_SEQUENCES = [
    "GYDPETGTWG",             # Chignolin (10 aa)
    "DAYAQWLKDGGPSSGRPPPS",   # Trp-cage (20 aa)
    "GEVEALKEKVSFLSALEEYT"    # Alpha Helix (20 aa)
]

# --- Gradio UI Creation Function ---
def create_demo():
    """Creates the Gradio interface application."""
    try: import py3Dmol
    except ImportError: logging.error("py3Dmol not found! 3D viewer will likely fail.")

    demo_app = BioEmuDemo() # Instantiate the class

    # Citation text
    citation_text_md = """
    ```bibtex
    @article {BioEmu2024,
        author = {Lewis, Sarah and Hempel, Tim and Jim{\'e}nez-Luna, Jos{\'e} and Gastegger, Michael and Xie, Yu and Foong, Andrew Y. K. and Satorras, Victor Garc{\'\i}a and Abdin, Osama and Veeling, Bastiaan S. and Zaporozhets, Iryna and Chen, Yaoyi and Yang, Soojung and Schneuing, Arne and Nigam, Jigyasa and Barbero, Federico and Stimper, Vincent and Campbell, Andrew and Yim, Jason and Lienen, Marten and Shi, Yu and Zheng, Shuxin and Schulz, Hannes and Munir, Usman and Clementi, Cecilia and No{\'e}, Frank},
        title = {Scalable emulation of protein equilibrium ensembles with generative deep learning},
        year = {2024},
        doi = {10.1101/2024.12.05.626885},
        journal = {bioRxiv}
    }
    ```
    """

    with gr.Blocks(title="BioEmu 3D Protein Demo", theme=gr.themes.Default()) as app:
        gr.Markdown("# BioEmu UI Wrapper - 3D Protein Structure Prediction Demo")
        gr.Markdown("""
        Generate 3D protein structures from amino acid sequences using the BioEmu model.\n
        Enter a sequence or load a sample. The resulting structure can be viewed interactively.
        """)

        output_dir_state = gr.State(None) # Stores the path to the output directory

        with gr.Row():
            with gr.Column(scale=3): # Input column
                input_text = gr.Textbox( label="Protein Sequence", placeholder="Enter amino acid sequence or FASTA format...", lines=5, elem_id="protein-sequence-input")
                with gr.Row():
                    num_samples_slider = gr.Slider(label="Number of Structures (BioEmu samples)", minimum=1, maximum=10, value=3, step=1, elem_id="num-samples-slider")
                with gr.Row():
                    submit_btn = gr.Button("Generate Structures", variant="primary", elem_id="generate-button")
                    sample_dropdown = gr.Dropdown(
                        label="Load Sample Sequence",
                        choices=[(f"Chignolin ({len(SAMPLE_SEQUENCES[0])} aa)", SAMPLE_SEQUENCES[0]),
                            (f"Trp-cage ({len(SAMPLE_SEQUENCES[1])} aa)", SAMPLE_SEQUENCES[1]),
                            (f"Alpha Helix ({len(SAMPLE_SEQUENCES[2])} aa)", SAMPLE_SEQUENCES[2])],
                        value=None, # No default selection
                        elem_id="sample-sequence-dropdown"
                    )
        with gr.Tabs(elem_id="output-tabs"):
            with gr.TabItem("Results Summary", elem_id="tab-summary"):
                output_text_box = gr.Textbox(label="Log & Summary", lines=10, interactive=False, show_copy_button=True, elem_id="summary-textbox")
            with gr.TabItem("3D Structure", elem_id="tab-structure"):
                # Using HTML component to embed the iframe for py3Dmol
                structure_viewer_html = gr.HTML(label="Structure Viewer (py3Dmol)", elem_id="structure-viewer-html")
            with gr.TabItem("Citation", elem_id="tab-citation"):
                gr.Markdown(citation_text_md, elem_id="citation-markdown")
            # Score plot is currently a placeholder and not generated by run_sample
            score_plot_output_plot = gr.Plot(visible=False, elem_id="score-plot") # Keep for future use

        # --- Download Area ---
        with gr.Row():
            download_btn = gr.Button("Download All Results (.zip)", elem_id="download-zip-button")
            download_file_output = gr.File(label="Download Zip File", interactive=False, elem_id="download-file-output")

        # --- Event Handlers ---
        # When a sample is chosen from the dropdown, update the input_text box
        sample_dropdown.change(fn=lambda choice_value: choice_value,
            inputs=sample_dropdown,
            outputs=input_text)

        # Wrapper function to handle progress for the BioEmu run
        def run_and_update_wrapper(seq, n_samples_from_slider, progress=gr.Progress(track_tqdm=True)):
            # Ensure n_samples is an integer for the bioemu_sample function
            n_samples = int(n_samples_from_slider)

            progress(0, desc="Initializing BioEmu...")
            logging.info(f"UI Trigger: Starting generation for sequence: {seq[:20]}..., N={n_samples}")
        
            # Visualise progress of BioEmu call
            progress(0.1, desc="Running BioEmu sampling...")
            summary, iframe_html, plot_out, out_dir = demo_app.run_sample(seq, n_samples)
            progress(0.9, desc="Processing results...")
        
            logging.info(f"run_sample completed. Dir: {out_dir}")
            progress(1, desc="Displaying results.")
            # The plot_out is currently always None from run_sample
            return summary, iframe_html, None, out_dir

        submit_btn.click(
            fn=run_and_update_wrapper,
            inputs=[input_text, num_samples_slider],
            outputs=[output_text_box, structure_viewer_html, score_plot_output_plot, output_dir_state]
        ).then(
            fn=lambda out_dir: gr.File(visible=False) if not out_dir else gr.File(visible=False), # Clear previous download link
            inputs=[output_dir_state],
            outputs=[download_file_output]
        )

        # When "Download Results" is clicked
        def download_results_interface(current_output_dir):
            if not current_output_dir:
                gr.Info("Please generate structures first before downloading.")
                return None # Return None if no directory to prevent errors
            zip_file_path = demo_app.download_structures(current_output_dir)
            if zip_file_path:
                gr.Info("Results zip file ready for download.")
                return gr.File(value=zip_file_path, visible=True) # Make file component visible with value
            else:
                # Warning/Error messages are handled within download_structures
                return gr.File(visible=False) # Hide if no file generated

        download_btn.click(
            fn=download_results_interface,
            inputs=[output_dir_state],
            outputs=[download_file_output]
        )

    logging.info("Gradio Blocks UI defined.")
    return app

# --- Main Execution Block ---
if __name__ == "__main__":
    log_format = '%(asctime)s - %(levelname)s - %(name)s - %(filename)s:%(lineno)d - %(message)s'
    logging.basicConfig(level=logging.INFO, format=log_format, force=True)
    logging.getLogger('PIL').setLevel(logging.WARNING) # Suppress less important PIL logs
    logging.getLogger('matplotlib').setLevel(logging.WARNING) # Suppress less important matplotlib logs
    logging.getLogger('huggingface_hub').setLevel(logging.INFO) # Can be WARNING in production

    logging.info("--- Starting BioEmu Gradio Demo Application ---")
    try:
        gradio_app = create_demo()
        logging.info("Gradio Demo UI created successfully.")
        # Launching the Gradio app
        logging.info("Launching Gradio interface on http://0.0.0.0:7860 (or the port you specify)")
        gradio_app.launch(
            server_name="0.0.0.0",
            server_port=7860,
            debug=True, # Enable Gradio debug mode for more detailed error messages in browser console
            # share=False # Set to True if you need a public link (requires internet)
        )
    except ImportError as ie:
        logging.exception(f"ImportError: {ie}. A required package might be missing. Check dependencies.")
        print(f"\n*** CRITICAL IMPORT ERROR: {ie} ***\nPlease ensure all required packages listed in requirements.txt are installed.\n")
    except NameError as ne:
        logging.exception(f"NameError: {ne}. This might be due to a typo or an undefined variable.")
        print(f"\n*** CRITICAL NAMERROR: {ne}. Please check the application code for typos or undefined names.\n")
    except TypeError as te:
        logging.exception(f"TypeError during setup or launch: {te}. Check function arguments or Gradio component configurations.")
        print(f"\n*** CRITICAL TYPEERROR: {te} ***\nThis often happens with incorrect arguments to functions or Gradio components.\n")
    except Exception as e:
        logging.exception("An unexpected critical error occurred during application setup or launch.")
        print(f"\n*** CRITICAL UNEXPECTED ERROR: {e} ***\n{traceback.format_exc()}\n")

    logging.info("--- BioEmu Gradio Demo Application Shutting Down ---")