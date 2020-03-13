import numpy as np
import matplotlib.pyplot as plt
from moviepy.editor import VideoClip
from moviepy.video.io.bindings import mplfig_to_npimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py
from matplotlib import rc
from matplotlib import ticker

# Plot settings
rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r'\usepackage{bm}']
rc('ps', usedistiller='xpdf')
rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
rc('axes', labelsize='xx-large')
rc('axes', grid='True')
rc('xtick', labelsize='xx-large')
rc('ytick', labelsize='xx-large')
plt.rcParams['axes.formatter.limits'] = (-2, 2)

# Create figures
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8.5, 8.5))
ax = axes.flatten()
divider1 = make_axes_locatable(ax[1])
divider2 = make_axes_locatable(ax[3])
cax1 = divider1.append_axes("right", "5%", pad="3%")
cax2 = divider2.append_axes("right", "5%", pad="3%")

# Input data
cbar1 = np.loadtxt('kp1e3_kb1e5_kn1e1_kg1e4_ka1e-3_ke1e3_kc1e3_M+005_Ns208_cbar.csv', delimiter=',', skiprows=1)
ctot1 = np.loadtxt('kp1e3_kb1e5_kn1e1_kg1e4_ka1e-3_ke1e3_kc1e3_M+005_Ns208_time.csv', usecols=9, delimiter=',', skiprows=1)
time1 = np.loadtxt('kp1e3_kb1e5_kn1e1_kg1e4_ka1e-3_ke1e3_kc1e3_M+005_Ns208_time.csv', usecols=0, delimiter=',', skiprows=1)
Pj1   =  h5py.File('kp1e3_kb1e5_kn1e1_kg1e4_ka1e-3_ke1e3_kc1e3_M+005_Ns208_prob_2D_t.jld2', 'r')

cbar2 = np.loadtxt('kp1e3_kb1e5_kn1e1_kg1e4_ka1e-3_ke1e-3_kc1e3_M+005_Ns208_cbar.csv', delimiter=',', skiprows=1)
ctot2 = np.loadtxt('kp1e3_kb1e5_kn1e1_kg1e4_ka1e-3_ke1e-3_kc1e3_M+005_Ns208_time.csv', usecols=9, delimiter=',', skiprows=1)
time2 = np.loadtxt('kp1e3_kb1e5_kn1e1_kg1e4_ka1e-3_ke1e-3_kc1e3_M+005_Ns208_time.csv', usecols=0, delimiter=',', skiprows=1)
Pj2   =  h5py.File('kp1e3_kb1e5_kn1e1_kg1e4_ka1e-3_ke1e-3_kc1e3_M+005_Ns208_prob_2D_t.jld2', 'r')

# Define Ns, diameter D (nm), i, and j
pref = 2.08
i = np.arange(1, cbar1.shape[1]+1)
Ns = np.around(pref*i**(2/3))
j = np.arange(0, Ns[-1]+1)
D = 0.25*(i/0.45)**(1/3)
I, J = np.meshgrid(i, j, indexing='ij')

# Calculate number density rho
P_j1 = Pj1['P_j']
rho1 = np.asarray([cbar1[t, :]/ctot1[t] for t in np.arange(0, time1.size)])
rho1[0, :] = 0.0

P_j2 = Pj2['P_j']
rho2 = np.asarray([cbar2[t, :]/ctot2[t] for t in np.arange(0, time2.size)])
rho2[0, :] = 0.0

# The solutions from different parameter sets do not have the same number of time points. Therefore, we fill
# the solutions with fewer time points with redundant data so that each solution has the same number of time points.
# This makes it simpler to make the movie.
threshold = -1.0
while time1.size < time2.size:
    diff = np.log(time2[:time1.size]) - np.log(time1)
    insert = np.where(diff <= threshold)[0]
    ind = np.ones(time1.size)
    ind[insert[0]] = 2
    rho1 = np.repeat(rho1, ind.astype(int), axis=0)
    P_j1 = np.repeat(P_j1, ind.astype(int), axis=0)
    time1 = np.insert(time1, insert[0], time1[insert[0]])

# Define function to pass to MoviePy
time1[0] = 0.0
time2[0] = 0.0


def make_frame(t):
    n = 2*t
    ax[0].clear()
    ax[1].clear()
    ax[2].clear()
    ax[3].clear()
    cax1.clear()
    cax2.clear()

    # Scheme I
    indices1 = np.asarray(np.where(rho1[int(n), :] >= 0.01))[0]

    ax[0].set_xlim(0.3, 1.4)
    if n == 0:
        ymax1 = 1.0
    else:
        ymax1 = rho1[int(n), indices1].max() + rho1[int(n), indices1].max()/10
    ax[0].set_ylim(0, ymax1)
    ax[0].set_ylabel(r'\textbf{Number Density}', labelpad=10.0)
    ax[0].plot(D, rho1[int(n), :])
    if indices1.size != 0:
        monomers1 = indices1 + 1
        x_coords1 = D[indices1]-0.015
        y_coords1 = rho1[int(n), indices1] + rho1[int(n), indices1]/30
        for index in np.arange(0, indices1.size):
            ax[0].text(x_coords1[index], y_coords1[index], monomers1[index].astype(str), fontsize='medium')

    plt.text(0.58, 1.05, 't = {0:.2e} s'.format(time1[int(n)]), transform=ax[0].transAxes, fontsize='x-large')

    ax[1].set_xlim(1, 68)
    ax[1].set_ylabel(r'$\bm{j}$')
    ax[1].set_ylim(0, Ns[-1])
    ax[1].set_title(r'$N_s = [2.08i^{2/3}]$', fontsize='x-large')
    cont1 = ax[1].contourf(I, J, P_j1[int(n), :, :].T, levels=4)
    ax[1].plot(i, Ns, color='black')

    colbar1 = plt.colorbar(cont1, cax=cax1)
    colbar1.locator = ticker.MaxNLocator(nbins=5)
    colbar1.update_ticks()
    colbar1.set_label(r'$\bm{P(j|i)}$', fontsize='x-large')
    plt.gcf().subplots_adjust(bottom=0.2, wspace=0.3)

    # Scheme IV
    indices2 = np.asarray(np.where(rho2[int(n), :] >= 0.01))[0]

    ax[2].set_xlim(0.3, 1.4)
    if n == 0:
        ymax2 = 1.0
    else:
        ymax2 = rho2[int(n), indices2].max() + rho2[int(n), indices2].max()/10
    ax[2].set_ylim(0, ymax2)
    ax[2].set_xlabel(r'$\bm{D}$ \textbf{(nm)}')
    ax[2].set_ylabel(r'\textbf{Number Density}', labelpad=10.0)
    ax[2].plot(D, rho2[int(n), :])
    if indices2.size != 0:
        monomers2 = indices2 + 1
        x_coords2 = D[indices2]-0.015
        y_coords2 = rho2[int(n), indices2] + rho2[int(n), indices2]/30
        for index in np.arange(0, indices2.size):
            ax[2].text(x_coords2[index], y_coords2[index], monomers2[index].astype(str), fontsize='medium')

    plt.text(0.58, 1.05, 't = {0:.2e} s'.format(time2[int(n)]), transform=ax[2].transAxes, fontsize='x-large')

    ax[3].set_xlim(1, 68)
    ax[3].set_xlabel(r'$\bm{i}$')
    ax[3].set_ylabel(r'$\bm{j}$')
    ax[3].set_ylim(0, Ns[-1])
    ax[3].set_title(r'$N_s = [2.08i^{2/3}]$', fontsize='x-large')
    cont2 = ax[3].contourf(I, J, P_j2[int(n), :, :].T, levels = 4)
    ax[3].plot(i, Ns, color='black')

    colbar2 = plt.colorbar(cont2, cax=cax2)
    colbar2.set_label(r'$\bm{P(j|i)}$', fontsize='x-large')
    colbar2.locator = ticker.MaxNLocator(nbins=5)
    colbar2.update_ticks()

    plt.gcf().subplots_adjust(bottom=0.2, wspace=0.3, hspace=0.3)

    return mplfig_to_npimage(fig)


animation = VideoClip(make_frame, duration=time1.size/2)
animation.write_videofile('movie_s1.mp4', fps=2)