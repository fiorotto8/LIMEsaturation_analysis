import ROOT
import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import pandas as pd
import os

import argparse

parser = argparse.ArgumentParser(description="Fit te saturation curves")
parser.add_argument("drift",help="put the DRIFT HV6 value",action="store",type=int)
parser.add_argument("-v", "--verbose", help="increase output verbosity",action="store_true")
args = parser.parse_args()

errIncrease=1.

def hist(data, x_name, channels=100, linecolor=4, linewidth=4,write=True):
    # Convert list directly to numpy array to avoid redundant loop
    array = np.array(data, dtype="d")
    # Create histogram
    hist = ROOT.TH1D(x_name, x_name, channels, 0.99*np.min(array), 1.01*np.max(array))
    # Use numpy vectorization to fill histogram
    for x in array:
        hist.Fill(x)
    # Set visual attributes and axis titles
    hist.SetLineColor(linecolor)
    hist.SetLineWidth(linewidth)
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle("Entries")
    # Set maximum digits on axes to manage display
    hist.GetYaxis().SetMaxDigits(3)
    hist.GetXaxis().SetMaxDigits(3)
    if write:
        hist.Write()
    return hist
def grapherr(x,y,ex,ey,x_string, y_string,name=None, color=4, markerstyle=22, markersize=2,write=True):
    plot = ROOT.TGraphErrors(len(x),  np.array(x  ,dtype="d")  ,   np.array(y  ,dtype="d") , np.array(   ex   ,dtype="d"),np.array( ey   ,dtype="d"))
    if name is None: plot.SetNameTitle(y_string+" vs "+x_string,y_string+" vs "+x_string)
    else: plot.SetNameTitle(name, name)
    plot.GetXaxis().SetTitle(x_string)
    plot.GetYaxis().SetTitle(y_string)
    plot.SetMarkerColor(color)#blue
    plot.SetMarkerStyle(markerstyle)
    plot.SetMarkerSize(markersize)
    if write==True: plot.Write()
    return plot
def split_string(s):
    try:
        before, after = s.split('/')
        return int(before), driftVtoEfield(int(after))
    except ValueError:
        print("Input string does not contain '/' or is not in the correct format.")
        return None, None
def draw_multigraph(graphs, title="MultiGraph", x_title="X-axis", y_title="Y-axis", names=None):
    """
    Draws a TMultiGraph on a canvas using a list of TGraphErrors,
    with automatically assigned marker styles, colors, and custom legend names.

    Parameters:
        graphs (list): List of TGraphErrors to be added to the TMultiGraph.
        title (str): Title of the TMultiGraph.
        x_title (str): X-axis title.
        y_title (str): Y-axis title.
        names (list): List of names for each graph (used in the legend).
    """
    # Create a canvas
    canvas = ROOT.TCanvas("canvas", title, 1000, 1000)
    canvas.SetLeftMargin(0.15)  # Adjust left margin to fit the Y-axis title

    # Create the TMultiGraph
    multi_graph = ROOT.TMultiGraph()
    multi_graph.SetTitle(f"{title};{x_title};{y_title}")

    # Predefined list of marker styles and colors for variety
    marker_styles = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29]  # ROOT marker styles
    colors = [2, 4, 6, 8, 9, 12, 28, 46, 30, 38]  # Different colors in ROOT

    # Loop over the graphs and assign styles and colors dynamically
    for i, graph in enumerate(graphs):
        # Cycle through marker styles and colors if we have more graphs than unique styles/colors
        style = marker_styles[i % len(marker_styles)]
        color = colors[i % len(colors)]

        # Apply marker style, color, and size for each graph
        graph.SetMarkerStyle(style)
        graph.SetMarkerColor(color)
        graph.SetMarkerSize(2)  # Adjust marker size as desired

        # Add the graph to the TMultiGraph
        multi_graph.Add(graph)

    # Draw the TMultiGraph
    multi_graph.Draw("AP")  # "A" for axis and "P" for points with error bars

    # Draw legend
    legend = ROOT.TLegend(0.15, 0.7, 0.25, 0.9)  # Adjust position as needed
    if names is None:
        names = [f"Graph {i+1}" for i in range(len(graphs))]  # Default names if not specified
    for i, graph in enumerate(graphs):
        legend.AddEntry(graph, names[i], "p")
    legend.Draw()

    # Update the canvas to display everything
    canvas.Update()
    canvas.Draw()
    canvas.SaveAs(f"{title}.png")
    canvas.Write()

    return canvas, multi_graph

def intTOgain(integral):
    conversionFactor=0.24 #e-/count
    lightFraction=0.07# gamma/e-
    omega=9.2E-4#solid angle fraction
    QE=0.7#QE @ 650nm
    n0=150#stima Fe55
    return (integral*conversionFactor)/(lightFraction*omega*QE*n0)
def sigmaTOdiff(sigma):
    pxsize = 155#um
    return float(sigma*pxsize)
def driftVtoEfield(V):
    d0_ring = 1.2#cm
    return int(V/d0_ring)
def soourcePostocm(x):
    return x*1.43

# Assuming args.drift has been defined earlier in the code
hv6drift = args.drift  # Set hv6drift to the value of args.drift

# Create the DataFrame from the CSV file
df = pd.read_csv("Position_Scan_Data_by_GEM_V_and_DRIFT_V_Configuration.csv")

# Filter columns based on hv6drift, keeping only those with the correct suffix after "/"
matching_columns = [col for col in df.columns if "/" in col and int(col.split("/")[1]) == hv6drift]
matching_columns.insert(0, "source_position")  # Add "source_position" to keep it in the DataFrame

# Create the filtered DataFrame and display it
filtered_df = df[matching_columns]
filtered_df = filtered_df.rename(columns={col: col.split("/")[0] for col in filtered_df.columns if "/" in col})

cuts_offset="(sc_rms>5) & (sc_tgausssigma>2.632) & (sc_width/sc_length>0.8) & (sc_integral<25000)& (sc_integral>1000) & (sc_tgausssigma<15)"
main= ROOT.TFile(f"out_{int(hv6drift)}Vcm.root","RECREATE")

graphs,names,vgems=[],[],[]
for vgem in filtered_df.columns[1:]:
    mean_integral, error_integral,mean_sigma, error_sigma,position = [], [], [], [], []
    main.mkdir(f"{vgem}")
    main.cd(f"{vgem}")
    vgems.append(vgem)
    for i, run in enumerate(filtered_df[vgem]):
        if args.verbose: print(vgem, run)
        if not pd.isna(run):
            pos=soourcePostocm(filtered_df["source_position"].iloc[i])
            file_name = f"reco_run{int(run)}_3D.root"
            if not os.path.exists("RECO/"+file_name):
                os.system(f"wget https://s3.cloud.infn.it/v1/AUTH_2ebf769785574195bde2ff418deac08a/cygno-analysis/RECO/Run4/{file_name}")
                os.system("mv "+file_name+" RECO/")
            try:
                events = uproot.open("RECO/"+file_name + ":Events")
            except:
                print("Failed to open (maybe empty)", file_name)
                continue

            cut_integral = events.arrays(["sc_integral"], f"{cuts_offset}")
            cut_sigma = events.arrays(["sc_tgausssigma"], f"{cuts_offset}")

            integralTEMP = [next(iter(d.values())) for d in ak.to_list(cut_integral)]
            sigmaTEMP = [next(iter(d.values())) for d in ak.to_list(cut_sigma)]
            
            integral=[item for sublist in integralTEMP for item in sublist]
            sigma=[item for sublist in sigmaTEMP for item in sublist]

            integral = [intTOgain(value) for value in integral]
            sigma = [sigmaTOdiff(value) for value in sigma]

            hist_integral = hist(integral, f"{vgem}_gain_{pos}",write=False)
            hist_sigma = hist(sigma, f"{vgem}_sigmaum_{pos}",write=False)
            
            range=[hist_integral.GetBinCenter(hist_integral.GetMaximumBin())-hist_integral.GetRMS(),hist_integral.GetBinCenter(hist_integral.GetMaximumBin())+hist_integral.GetRMS()]
            
            #integral
            gaus = ROOT.TF1("gaus", "gaus(0)", range[0],range[1])
            gaus.SetParameters(hist_integral.GetMaximum(), hist_integral.GetMean(), hist_integral.GetRMS())
            hist_integral.Fit(gaus, "RQ")
            mean_integral.append(gaus.GetParameter(1))        
            error_integral.append(errIncrease*gaus.GetParError(1))
            hist_integral.Write()
            #tgaussigma
            gaus = ROOT.TF1("gaus", "gaus(0)", 100,2000)
            gaus.SetParameters(hist_sigma.GetMaximum(), hist_sigma.GetMean(), hist_sigma.GetRMS())
            hist_sigma.Fit(gaus, "RQ")
            mean_sigma.append(gaus.GetParameter(1))        
            error_sigma.append(errIncrease*gaus.GetParError(1))
            hist_sigma.Write()
            
            position.append(pos)

    graphs.append(grapherr(mean_sigma, mean_integral, error_sigma, error_integral, "sc_tgausssigma", "sc_integral", f"{vgem}_integral_vs_sigma",write=False))
    #graphs.append(grapherr(position, mean_integral, np.zeros(len(position)), error_integral, "sc_tgausssigma", "sc_integral", f"{vgem}_integral_vs_sigma",write=False))
    names.append(vgem)

if args.verbose: print(vgems)

main.cd()

rows=[]

for i,graph in enumerate(graphs): 
    #Gtot = ROOT.TF1("Gtot", f"[0]^3 * x^3 / (x^3 + ([1]/{vgems[i]}) * [0]^2 * ([0] - 1))", 0, 2000)
    #GtotAtt = ROOT.TF1("GtotAtt", f"([0]^3 * x^3 / (x^3 + ([1]/{vgems[i]}) * [0]^2 * ([0] - 1))) * (exp(-[2] * (x*x - [3])))", 0, 2000)

    Gtot = ROOT.TF1("Gtot", f"(exp(3*[0]*{vgems[i]}) * x^3) / (x^3 + ([1]/{vgems[i]}) * exp(2*[0]*{vgems[i]}) * (exp([0]*{vgems[i]}) - 1))", 0, 2000)
    GtotAtt = ROOT.TF1("GtotAtt", f"(exp(3*[0]*{vgems[i]}) * x^3) / (x^3 + ([1]/{vgems[i]}) * exp(2*[0]*{vgems[i]}) * (exp([0]*{vgems[i]}) - 1))* (exp(-[2] * (x*x - [3])))", 0, 2000)

    Gtot.SetParameters(50, 1.3E5)
    Gtot.SetParameters(0.009, 1.3E5)
    Gtot.FixParameter(1,1.3E5)
    graph.Fit(Gtot, "RQ")
    Gtot.ReleaseParameter(1)
    graph.Fit(Gtot, "RQ")

    GtotAtt.SetParameters(Gtot.GetParameter(0), Gtot.GetParameter(1), 4E-7, 1E6)
    GtotAtt.FixParameter(1,Gtot.GetParameter(1))
    GtotAtt.FixParameter(0,Gtot.GetParameter(0))
    #GtotAtt.FixParameter(2,1E-7)
    #GtotAtt.FixParameter(3,900)
    graph.Fit(GtotAtt, "RQ")
    
    GtotAtt.ReleaseParameter(1)
    GtotAtt.ReleaseParameter(0)
    graph.Fit(GtotAtt, "RQ")

    # Define the row you want to append as a string
    rows.append( f"{driftVtoEfield(hv6drift)};{vgems[i]};" \
                f"{GtotAtt.GetParameter(0):.2e};{GtotAtt.GetParError(0):.2e};" \
                f"{GtotAtt.GetParameter(1):.2e};{GtotAtt.GetParError(1):.2e};" \
                f"{GtotAtt.GetParameter(2):.2e};{GtotAtt.GetParError(2):.2e};" \
                f"{GtotAtt.GetParameter(3):.2e};{GtotAtt.GetParError(3):.2e};" \
                f"{GtotAtt.GetChisquare()/GtotAtt.GetNDF():.2f}\n" )
                
    # Open the file in append mode and write the new row
    with open("parameters.txt", "a") as file: 
        file.write(rows[i])
    
    graph.Write()

    if args.verbose: print("First fit")
    if args.verbose: print(graph.GetName(),"p0:",f"{Gtot.GetParameter(0):.2e}","p1:",f"{Gtot.GetParameter(1):.2e}", "chi2/ndf",f"{(Gtot.GetChisquare()/Gtot.GetNDF()):.2f}")
    if args.verbose: print("Second fit")
    if args.verbose: print(graph.GetName(),"p0:",f"{GtotAtt.GetParameter(0):.2e}","p1:",f"{GtotAtt.GetParameter(1):.2e}","p2:",f"{GtotAtt.GetParameter(2):.2e}","p3:",f"{GtotAtt.GetParameter(3):.2e}", "chi2/ndf",f"{(GtotAtt.GetChisquare()/GtotAtt.GetNDF()):.2f}")


draw_multigraph(graphs, title=f"Gain vs Sigma VDRIFT {driftVtoEfield(hv6drift)} Vcm", x_title="Sigma (#mum)", y_title="Gain", names=names)
