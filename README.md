# proteosim_JKa
## What does this project do?

This project simulates a simplified bottom-up proteomics workflow. It covers several key steps commonly used in protein analysis:

- **Protein digestion** into peptide fragments
- **Liquid chromatography (LC)**-based peptide separation using predicted retention times
- **Mass spectrometry (MS)** characterization by calculating peptide masses and m/z values
- **MS/MS fragmentation** by generating b- and y-ion series and calculating their corresponding masses and m/z values
- **Data visualization** through retention time distributions and simulated mass spectra
- **Export of results** as TSV files for further analysis

The project can be used as an educational tool to understand the computational principles behind a proteomics pipeline.

## What are the requirements or installation steps?

### Requirements

- Python 3.10 or newer
- NumPy
- Pandas
- Matplotlib
- Pyteomics

### Installation

Clone the repository:

git bash
git clone <repository-url>
cd <repository-name>
