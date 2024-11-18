import numpy as np
import ROOT

# Function to read data from parameters.txt file
def read_parameters(file_path):
    data = np.genfromtxt(file_path, delimiter=';', skip_header=1)
    drift_fields = data[:, 0]
    vgem = data[:, 1]
    p0 = data[:, 2]
    ep0 = data[:, 3]
    p1 = data[:, 4]
    ep1 = data[:, 5]
    p2 = data[:, 6]
    ep2 = data[:, 7]
    p3 = data[:, 8]
    ep3 = data[:, 9]
    chi2ndf = data[:, 10]
    return drift_fields, vgem, p0, ep0, p1, ep1, p2, ep2, p3, ep3,chi2ndf
# Your grapherr function as given
def grapherr(x, y, ex, ey, x_string, y_string, name=None, color=4, markerstyle=22, markersize=2, write=True):
    plot = ROOT.TGraphErrors(len(x), np.array(x, dtype="d"), np.array(y, dtype="d"), np.array(ex, dtype="d"), np.array(ey, dtype="d"))
    if name is None:
        plot.SetNameTitle(y_string + " vs " + x_string, y_string + " vs " + x_string)
    else:
        plot.SetNameTitle(name, name)
    plot.GetXaxis().SetTitle(x_string)
    plot.GetYaxis().SetTitle(y_string)
    plot.SetMarkerColor(color)
    plot.SetMarkerStyle(markerstyle)
    plot.SetMarkerSize(markersize)
    if write:
        plot.Write()
    return plot
# Draw multiple graphs on a single canvas
def draw_multigraph(graphs, title="MultiGraph", x_title="X-axis", y_title="Y-axis", names=None):
    canvas = ROOT.TCanvas("canvas", title, 1000, 1000)
    canvas.SetLeftMargin(0.15)

    multi_graph = ROOT.TMultiGraph()
    multi_graph.SetTitle(f"{title};{x_title};{y_title}")

    marker_styles = [20, 21, 22, 23, 33,34]
    colors = [2, 4, 6, 8, 9, 12, 28, 46, 30, 38]

    for i, graph in enumerate(graphs):
        style = marker_styles[i % len(marker_styles)]
        color = colors[i % len(colors)]
        graph.SetMarkerStyle(style)
        graph.SetMarkerColor(color)
        graph.SetMarkerSize(2)
        multi_graph.Add(graph)

    multi_graph.Draw("AP")

    legend = ROOT.TLegend(0.15, 0.7, 0.4, 0.9)
    if names is None:
        names = [f"Graph {i+1}" for i in range(len(graphs))]
    for i, graph in enumerate(graphs):
        legend.AddEntry(graph, names[i], "p")
    legend.Draw()

    canvas.Update()
    canvas.Draw()
    canvas.SaveAs(f"{title}.png")
    canvas.Write()

    return canvas, multi_graph


# Main script to read the data and plot each parameter
file_path = 'parameters.txt'  # Adjust the path if necessary
drift_fields, vgem, p0, ep0, p1, ep1, p2, ep2, p3, ep3, chi2ndf = read_parameters(file_path)

root_file = ROOT.TFile("output_parameter_plots.root", "RECREATE")

# Get unique values of Drift_field and VGEM
unique_drift_fields = np.unique(drift_fields)
unique_vgem = np.unique(vgem)

# Parameters and their errors
parameters = [
    (p0, ep0, "p0"), (p1, ep1, "p1"), 
    (p2, ep2, "p2"), (p3, ep3, "p3"),
    (chi2ndf, np.zeros_like(chi2ndf), "chi2ndf")
]

# Plot each parameter as a function of VGEM with multiple Drift Fields on each canvas
for param_data, param_error, param_name in parameters:
    graphs = []
    names = []
    
    for drift_field in unique_drift_fields:
        indices = np.where(drift_fields == drift_field)
        vgem_filtered = vgem[indices]
        param_filtered = param_data[indices]
        param_error_filtered = param_error[indices]
        
        # Create graph for the current drift field and parameter
        graph = grapherr(vgem_filtered, param_filtered, np.zeros_like(vgem_filtered), param_error_filtered, "VGEM (V)", param_name, color=4)
        graphs.append(graph)
        names.append(f"{param_name} (Drift Field = {int(drift_field)} V/cm)")
    
    # Draw all drift field graphs for the current parameter on one canvas
    title = f"{param_name} vs VGEM for Different Drift Fields"
    draw_multigraph(graphs, title=title, x_title="VGEM (V)", y_title=param_name, names=names)

# Plot each parameter as a function of Drift Field with multiple VGEMs on each canvas
for param_data, param_error, param_name in parameters:
    graphs = []
    names = []
    
    for vgem_value in unique_vgem:
        indices = np.where(vgem == vgem_value)
        drift_filtered = drift_fields[indices]
        param_filtered = param_data[indices]
        param_error_filtered = param_error[indices]
        
        # Create graph for the current VGEM value and parameter
        graph = grapherr(drift_filtered, param_filtered, np.zeros_like(drift_filtered), param_error_filtered, "Drift Field (V/cm)", param_name, color=4)
        graphs.append(graph)
        names.append(f"{param_name} (VGEM = {int(vgem_value)} V)")
    
    # Draw all VGEM graphs for the current parameter on one canvas
    title = f"{param_name} vs Drift Field for Different VGEMs"
    draw_multigraph(graphs, title=title, x_title="Drift Field (V/cm)", y_title=param_name, names=names)

root_file.Close()