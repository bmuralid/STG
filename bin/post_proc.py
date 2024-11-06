from plotsettings import *

def plot_sem_rans_comp(fpath='./'):
    uprime = np.loadtxt(f'{fpath}/uprime.dat')
    vprime = np.loadtxt(f'{fpath}/vprime.dat')
    wprime = np.loadtxt(f'{fpath}/wprime.dat')
    
    dat_cord = np.loadtxt('stg_interface_01.dat', skiprows=1)
    
    y = dat_cord[:, 1]
    
    uu = uprime * uprime /1000.0
    vv = vprime * vprime /1000.0
    ww = wprime * wprime /1000.0
    uv = uprime * vprime /1000.0
    
    tke = 0.5 * (uu + vv + ww)
    
    dat = np.loadtxt('wall_prof.dat', skiprows=1)
    
    tke_avg = []
    urms_avg = []
    vrms_avg = []
    wrms_avg = []
    uvrms_avg = []
    yloc = []
    for i in range(np.shape(dat)[0]):
        #indx = np.where(np.abs(y[:]-dat[i, 0]) < 0.00001*dat[i,0]) 
        indx = np.where(np.abs(y[:]-dat[i, 0]) < 0.001*dat[i,0]) 
        if len(indx[0]) > 0:
            tke_avg.append(np.sum(tke[indx])/len(indx[0]))
            urms_avg.append(np.sum(uu[indx])/len(indx[0]))
            vrms_avg.append(np.sum(vv[indx])/len(indx[0]))
            wrms_avg.append(np.sum(ww[indx])/len(indx[0]))
            uvrms_avg.append(np.sum(uv[indx])/len(indx[0]))
            yloc.append(dat[i,0])
    
    tke_avg = np.array(tke_avg)
    urms_avg = np.array(urms_avg)
    vrms_avg = np.array(vrms_avg)
    wrms_avg = np.array(wrms_avg)
    uvrms_avg = np.array(uvrms_avg)
    yloc = np.array(yloc)
    
    yy = dat[:, 0]
    deltay = dat[:, 1]
    u = dat[:, 2]
    tke_rans = dat[:, 3]
    epsilon = dat[:, 4]
    
    delta0 = 1.0
    
    plt.figure()
    plt.semilogx(yloc/delta0,tke_avg,'-b',label='SEM-TKE')
    plt.semilogx(yy/delta0,tke_rans,'-k', label='RANS-TKE')
    plt.legend(loc='best')
    plt.xlabel(r'$y [m]$', fontsize=24)
    plt.ylabel('TKE [m^s/s^2]')
    #plt.savefig('tke_comp_sem.png', dpi=200, bbox_inches='tight')
    
    R11 = 4.0/9.0 * 2.0 * tke_rans
    R22 = 2.0/9.0 * 2.0 * tke_rans
    R33  = 1.0/3.0 * 2.0 * tke_rans
    R12 = -0.3*tke_rans
    
    plt.figure()
    plt.semilogx(yloc/delta0, urms_avg, '-k', label='<uu>')
    plt.semilogx(yy/delta0, R11, '--k')
    
    #plt.figure()
    plt.semilogx(yloc/delta0, vrms_avg, '-b', label='<vv>')
    plt.semilogx(yy/delta0, R22, '--b')
    
    #plt.figure()
    plt.semilogx(yloc/delta0, wrms_avg, '-r', label='<ww>')
    plt.semilogx(yy/delta0, R33,'--r')
    
    #plt.figure()
    plt.semilogx(yloc/delta0, uvrms_avg, '-g', label='<uv>')
    plt.semilogx(yy/delta0, R12, '--g')
    plt.xlabel(r'$y [m]$', fontsize=24)
    plt.ylabel(r'$<u_iu_j> [m^s/s^2]$', fontsize=24)
    ##plt.savefig('rij_comp_sem.png', dpi=200, bbox_inches='tight')
    
    plt.show()


if __name__ == '__main__':
    plot_sem_rans_comp('./')
