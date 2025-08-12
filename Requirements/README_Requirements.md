# 📦 ASO Pipeline - Requirements

This folder contains all the necessary information to set up the Python environment for the **ASO Design and Exon Skipping Simulation Pipeline**.

## 🧪 Required Packages

All dependencies are listed in the `requirements.txt` file.

### 📋 Installation Instructions

To install all required packages, run the following command from the root of the repository:

```bash
pip install -r Requirements/requirements.txt
This will ensure that the pipeline runs smoothly with the correct versions of:
```

## 📚 Package Descriptions

Below is an overview of the Python packages used in this project and their roles.

### 1. BioPython
A collection of tools for biological computation, including sequence parsing, manipulation, and file format conversion (e.g., FASTA, GenBank, GTF). Essential for handling genomic and proteomic data in the pipeline.

### 2. pandas
A fast and flexible data analysis library, providing powerful DataFrame structures for storing, filtering, and processing tabular results such as exon annotations, ASO candidates, and simulation outputs.

### 3. numpy
A fundamental package for numerical and scientific computing in Python. It supports efficient array operations, which are critical for sequence indexing, coordinate calculations, and statistical analysis.

### 4. matplotlib
A widely used plotting library for creating static, animated, and interactive visualizations. In this pipeline, it is used for visualizing exon structures, skipping patterns, and analytical results.

### 5. scikit-learn
A machine learning library offering tools for classification, regression, and clustering. While optional in the current pipeline, it can be used for advanced ASO scoring, feature ranking, and predictive modeling.

### 6. requests
A simple yet powerful HTTP library for sending requests to APIs and downloading online resources (e.g., retrieving Ensembl or UniProt data).

### 7. tqdm
A lightweight progress bar library that provides real-time feedback during long-running tasks, such as parsing large GTF files or iterating through all possible exon skip combinations.

### 8. beautifulsoup4
A library for parsing HTML and XML documents. In this project, it is primarily used to extract protein domain annotations and metadata from UniProt XML files.

### 9. lxml
A high-performance XML and HTML parsing library used in conjunction with BeautifulSoup for faster and more efficient data extraction.

### 10. seaborn
A statistical data visualization library based on matplotlib, offering more attractive and informative graphics. It is useful for heatmaps, distribution plots, and comparative analyses of ASO candidates.



✅ If you're using a virtual environment (recommended), activate it before running the install command.

##💡 Tip
### To create and activate a virtual environment:

python -m venv aso_env
source aso_env/bin/activate  # On Windows: aso_env\Scripts\activate
pip install -r Requirements/requirements.txt
