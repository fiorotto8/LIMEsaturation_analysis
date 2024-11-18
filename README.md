# Anlayse saturation data from LIME campaign

- `ListRun.csv` contain the segmetn of the run log wiht the runs used
- `Position_Scan_Data_by_GEM_V_and_DRIFT_V_Configuration.csv` is a file fgoruping the runs for VGEM and DRIFT_V

## `anal.py` 

This script analyzes GEM detector data by fitting saturation curves. It filters and processes event data, applies cuts, and generates histograms and graphs of various metrics such as gain and sigma values.

### Usage

To execute the script, provide the required drift voltage value (DRIFT HV6) as an argument and an optional verbosity flag for detailed output.

```bash
python script.py <drift_voltage> [-v]
```

- `<drift_voltage>`: Integer, DRIFT HV6 value.
- `-v, --verbose`: Optional flag to increase output verbosity.

### Workflow

The script uses two main fitting functions applied to the gain versus sigma data for each GEM voltage configuration. These functions aim to model the gain as a function of sigma, which relates to the behavior of the GEM detector under different conditions.

1. **General Gain Function (`Gtot`)**:
   \[
   G_{\text{tot}}(x) = \frac{\exp(3a \cdot V_{\text{GEM}}) \cdot x^3}{x^3 + \left(\frac{b}{V_{\text{GEM}}}\right) \cdot \exp(2a \cdot V_{\text{GEM}}) \cdot (\exp(a \cdot V_{\text{GEM}}) - 1)}
   \]


2. **Gain Function with Attenuation (`GtotAtt`)**:
   \[
   G_{\text{totAtt}}(x) = G_{\text{tot}}(x) \cdot \exp\left(-c \cdot (x^2 - d)\right)
   \]


These functions are fitted to the gain and sigma data for each GEM voltage configuration to capture the behavior of gain under varying detector conditions. The fitted parameters \( a \), \( b \), \( c \), and \( d \) are saved to `parameters.txt` for further analysis.

### Output

- Root files and graphs in `.png` format are saved based on the input drift voltage.
- `parameters.txt` file stores the fit parameters for further analysis.

## `results.py` 

This script reads and analyzes GEM detector parameters from the `parameters.txt` file, creating plots to visualize the relationship between different parameters (such as gain and sigma) as functions of GEM voltage (VGEM) and drift field. The output includes both ROOT and PNG files with parameter plots.

### Requirements

- Python 3
- ROOT (CERN)
- NumPy

### Usage

Ensure the `parameters.txt` file is in the same directory or specify the path in `file_path`

### Plot Types

The script produces two sets of plots for each parameter (p0, p1, p2, p3, and chi2ndf):
1. **Parameter vs VGEM**: Each drift field configuration is plotted on the same canvas for a given parameter as a function of VGEM.
2. **Parameter vs Drift Field**: Each VGEM configuration is plotted on the same canvas for a given parameter as a function of drift field.

### Output

- **ROOT File**: `output_parameter_plots.root` contains all the plots generated.
- **PNG Images**: Each plot is also saved as a PNG file for easy viewing.
