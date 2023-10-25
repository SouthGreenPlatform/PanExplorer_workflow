import numpy as np
import pandas as pd
import seaborn as sns
import plotly.express as px
import xarray as xr

import sys, getopt

np.bool = np.bool_


def main(argv):
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print ('Heatmap.py -i <inputfile> -o <outputfile>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('Heatmap.py -i <inputfile> -o <outputfile>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
   print ('Input file is "', inputfile)
   print ('Output file is "', outputfile)

   data = pd.read_csv(inputfile, sep='\t')

   #data.drop(columns=data.columns[0], axis=1,  inplace=True)
   
   data2 = pd.read_csv(inputfile, sep='\t', header=None)
   data2.drop(index=data2.index[0], axis=0, inplace=True)
   data3 = data2[data2.columns[1:]]
   data3_transposed = data3.T

   myList = list(data.columns)
   del myList[0]

   myClusters = list(data['ClutserID'])

   fig = px.imshow(data3_transposed,
                x=myClusters,
                y=myList,
                labels=dict(y="Strains", x="Clusters"),
                height = 900,width = 900,
                color_continuous_scale=["lightgrey", "red","green"]
               )

   fig.update_traces(hovertemplate="<br>".join(["Cluster: %{x}","Strain: %{y}"]))
   fig.update_coloraxes(showscale=False)
   fig.update_layout(xaxis_scaleanchor="x")
   fig.write_html(outputfile)

if __name__ == "__main__":
   main(sys.argv[1:])
