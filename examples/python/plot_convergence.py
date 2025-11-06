import numpy as np
import matplotlib.pyplot as plt


# Configure 
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

# Optionally set font to Computer Modern to avoid common missing
# font errors
params = {
  'axes.labelsize': 16,
  'legend.fontsize': 14,
  'xtick.labelsize': 16,
  'ytick.labelsize': 16,
  'text.usetex': True}
plt.rcParams.update(params)

# Latex math
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath}\usepackage{amsmath}\usepackage{amssymb}']
#plt.rcParams['font.family'] = 'sans-serif'
#plt.rcParams['font.sans-serif'] = 'courier'
plt.rcParams['font.size'] = 16
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['lines.linewidth'] = 2
plt.rcParams['lines.color'] = 'r'

# Make sure everything is within the frame
plt.rcParams.update({'figure.autolayout': True})

# Set marker size
markerSize = 7.0 #11.0
mew = 2.0

# bar chart settings
lalpha    = 0.9
rev_alpha = 0.9/1.5

# These are the "Tableau 20" colors as RGB.    
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (25, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229),
             (255,223,0)]    
  
# Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)
    
plt.figure()
fig, ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
   
# Define colors
poscolor = tableau20[2]
posbandcolor = tableau20[3]

velcolor = tableau20[16]
velbandcolor = tableau20[17]

acccolor = tableau20[0]
accbandcolor = tableau20[1]

alpha_level = 0.3

# COMPLETE data
complete_terms = np.array([1,4,10,20,35,56,84,120])
complete_mean_error = np.array([0.017154474598167282,
                                0.0003468756341509456,
                                1.0976874260729254e-05,
                                5.615464522786541e-07,
                                2.3688839313538083e-08,
                                1.1163844549344318e-09 ,
                                7.743791886305301e-11,
                                7.593252074049514e-12 ])
complete_variance_error = np.array([0.0012732482140375946,
                               0.00010034804241925243,
                               2.6307903349504656e-06,
                               1.3866247660018434e-07,
                               9.498226899007596e-09,
                               5.539090935511577e-10,
                               3.8854368380976567e-11,
                               3.74712269662585e-12])

# 216 1.1436804032072324e-09 5.791400056404112e-10

# TENSOR data
tensor_terms = np.array([1,8,27,64,125,216])
tensor_mean_error = np.array([ 0.017154474598167282,
                            0.00034806571941763044,
                            1.1175375459382694e-05,
                            5.746518374944619e-07,
                            2.4298843790143135e-08,
                            1.1436804032072324e-09])
tensor_variance_error = np.array([0.0012732482140375946,
                                  9.189659752193513e-05,
                                  2.4208285005551143e-06,
                                  1.4270097438540776e-07,
                                  9.830884863487853e-09,
                                  5.791400056404112e-10])

plt.loglog(complete_terms, complete_mean_error, '-o'    , label = 'expectation~:~complete~basis', color=tableau20[14], mew=mew, ms=markerSize, mec='black')
plt.loglog(tensor_terms, tensor_mean_error, '-s'        , label = 'expectation~:~tensor~basis'  , color=tableau20[14], mew=mew, ms=markerSize, mec='black')
plt.loglog(complete_terms, complete_variance_error, '-o', label = 'variance~~~~~:~complete~basis', color=tableau20[16], mew=mew, ms=markerSize, mec='black')
plt.loglog(tensor_terms, tensor_variance_error, '-s'     ,label = 'variance~~~~~:~tensor~basis', color=tableau20[16], mew=mew, ms=markerSize, mec='black')
plt.legend(loc='upper right', frameon=False, ncol=1, handletextpad=0.3)
plt.xlabel('cardinality of orthonormal basis set')
plt.ylabel('norm of absolute error')
plt.ylim(top=1.0e-1,bottom=1.0e-12)
plt.savefig('smd-convergence-complete-tensor-loglog.pdf', bbox_inches='tight', pad_inches = 0.035)
plt.close()

plt.figure()
fig, ax = plt.subplots()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.semilogy()

plt.semilogy(complete_terms, complete_mean_error    , '-o' , label = 'expectation~:~complete~basis', color=tableau20[14], mew=mew, ms=markerSize/1.5, mec='black',alpha=0.750)
plt.semilogy(tensor_terms, tensor_mean_error    , '-s'     , label = 'expectation~:~tensor~basis'  , color=tableau20[14], mew=mew, ms=markerSize/1.5, mec='black',alpha=0.750)
plt.semilogy(complete_terms, complete_variance_error, '-o' , label = 'variance~~~~~:~complete~basis'   , color=tableau20[16], mew=mew, ms=markerSize/1.5, mec='black',alpha=0.750)
plt.semilogy(tensor_terms, tensor_variance_error, '-s'     , label = 'variance~~~~~:~tensor~basis'     , color=tableau20[16], mew=mew, ms=markerSize/1.5, mec='black',alpha=0.750)
#plt.legend(loc='upper left', frameon=False, ncol=2, handletextpad=0.3,  bbox_to_anchor=(0.0, 1.07))
plt.legend(loc='upper right', frameon=False, ncol=1, handletextpad=0.3)
plt.xlabel('cardinality of orthonormal basis set')
plt.ylabel('norm of absolute error')
plt.ylim(top=1.0,bottom=1.0e-12)
plt.xticks([0,30,60,90,120,150,180,210])
plt.savefig('smd-convergence-complete-tensor-semilogy.pdf', bbox_inches='tight', pad_inches = 0.035)
