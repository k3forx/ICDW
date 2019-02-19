import configparser
import numpy as np
########################################
########### read parameter #############
########################################
config = configparser.ConfigParser()
config.read('pHMC_PARAM')

NSITE = config.getint('param','NSITE')
dim = config.getint('param','dim')

mchi = config.getfloat('param','mchi')
gg = config.getfloat('param','gg')
lam = config.getfloat('param','lam')
nu = config.getfloat('param','nu')

NTHERM = config.getint('param','NTHERM')
NSKIP = config.getint('param','NSKIP')
NSAMPLE = config.getint('param','NSAMPLE')

nmd = config.getint('param','nmd')
iseed = config.getint('param','seed')
########################################################
########### set up of unit vector and periodic vector
########################################################
mu = np.zeros((dim,dim),dtype='int64')
for i in range(dim):
    mu[i,i] = 1

pn = np.zeros((NSITE),dtype='int64')
for i in range(NSITE): # i = 0...NSITE-1
    pn[i] = i
pn = np.insert(pn,NSITE,[0,1])
pn = np.insert(pn,0,[NSITE-2,NSITE-1])

################################################
# computation of hamiltonian
################################################
def hamil (sx,sp,vx,vp,mchi,lam,nu,gg):
    h = 0.0

    for i in range(NSITE):
        for j in range(NSITE):
            sum_v = sum_vec(vp,i,j)
            div_s = div_sca(sx,i,j)

            h = h + 0.5*sp[i,j]**2 - 0.5*sx[i,j]*div_s + 0.5*sum_v + lam*(sx[i,j]**2 - nu**2)**2

            sum_v = sum_vec(vx,i,j)
            div_v = div_vec(vx,i,j)

            h = h + 0.5*sum_v + 0.5*(div_v - gg*sx[i,j])**2/mchi**2

    return h


#######################################################
######## leapfrog molecular dynamics
#######################################################
def leapfrog_xpx_md (tau,nmd,sx,sp,vx,vp,mchi,lam,nu,gg):
    dt = tau/nmd
    dt2 = dt/2.0

    sx,vx = update_x(dt2,sx,sp,vx,vp,mchi,lam,nu,gg)

    for i in range(nmd-1):
        sp,vp = update_p(dt,sx,sp,vx,vp,mchi,lam,nu,gg)
        sx,vx = update_x(dt,sx,sp,vx,vp,mchi,lam,nu,gg)

    sp,vp = update_p(dt,sx,sp,vx,vp,mchi,lam,nu,gg)
    sx,vx = update_x(dt2,sx,sp,vx,vp,mchi,lam,nu,gg)

    return sx,sp,vx,vp

def update_x (dt,sx,sp,vx,vp,mchi,lam,nu,gg):
    return sx+sp*dt,vx+vp*dt


def update_p (dt,sx,sp,vx,vp,mchi,lam,nu,gg):
    sf = get_scalar_force(sx,vx,mchi,lam,nu,gg)
    vf = get_vector_force(sx,vx,mchi,lam,nu,gg)

    return sp+sf*dt,vp+vf*dt


######################################################
###### calculation of force
######################################################
def get_scalar_force (sx,vx,mchi,lam,nu,gg):
    sf = np.zeros((NSITE,NSITE),dtype='float64')

    for i in range(NSITE):
        for j in range(NSITE):
            div_s = div_sca(sx,i,j)
            div_v = div_vec(vx,i,j)
            sf[i,j] = div_s - 4.0*lam*sx[i,j]*(sx[i,j]**2 - nu**2) \
                        + gg*(div_v - gg*sx[i,j])/mchi**2
    return sf


def get_vector_force (sx,vx,mchi,lam,nu,gg):
    vf = np.zeros((NSITE,NSITE,dim),dtype='float64')

    for k in range(dim):
        for j in range(NSITE):
            for i in range(NSITE):
                pni = pn[i + mu[k,0] + 2]
                pnj = pn[j + mu[k,1] + 2]
                div_v = div_vec(vx,pni,pnj)

                vf[i,j,k] = - vx[i,j,k] + 0.5*(div_v - gg*sx[pni,pnj])/mchi**2

                pni = pn[i - mu[k,0] + 2]
                pnj = pn[j - mu[k,1] + 2]
                div_v = div_vec(vx,pni,pnj)

                vf[i,j,k] = + vf[i,j,k] - 0.5*(div_v - gg*sx[pni,pnj])/mchi**2

    return vf


######################################################
####### derivative part
######################################################
def div_sca (s,i,j):
    div_s = 0.0
    for k in range(dim):
        div_s += s[pn[i+mu[k,0]+2],pn[j+mu[k,1]+2]] \
                    - 2.0*s[i,j] + s[pn[i-mu[k,0]+2],pn[j-mu[k,1]+2]]
    return div_s


def div_vec (v,i,j) :
    div_v = 0.0
    for k in range(dim) :
        div_v += v[pn[i+mu[k,0]+2],pn[j+mu[k,1]+2],k] - v[pn[i-mu[k,0]+2],pn[j-mu[k,1]+2],k]
    return 0.5*div_v


def sum_vec (v,i,j):
    sum_v = 0.0
    for k in range(dim):
        sum_v += v[i,j,k]**2
    return sum_v


################################################
################################################
################ main program ##################
################################################
################################################
np.random.seed(seed=iseed)

sx0 = np.random.normal(0,1,(NSITE,NSITE))
vx0 = np.random.normal(0,1,(NSITE,NSITE,dim))

tau = 1.0
itry = 0
iacc = 0

for istep in range(NTHERM + NSKIP*NSAMPLE):
    sp0 = np.random.normal(0,1,(NSITE,NSITE))
    vp0 = np.random.normal(0,1,(NSITE,NSITE,dim))

    h0 = hamil(sx0,sp0,vx0,vp0,mchi,lam,nu,gg)
    # print(h0)

    sx1 = sx0
    sp1 = sp0
    vx1 = vx0
    vp1 = vp0

    sx1,sp1,vx1,vp1 = leapfrog_xpx_md(tau,nmd,sx1,sp1,vx1,vp1,mchi,lam,nu,gg)

    sp1 = -sp1
    vp1 = -vp1

# _CHECK_REVERSE_
    sx2 = sx1
    sp2 = sp1
    vx2 = vx1
    vp2 = vp1

    sx2,sp2,vx2,vp2 = leapfrog_xpx_md(tau,nmd,sx2,sp2,vx2,vp2,mchi,lam,nu,gg)
    h2 = hamil(sx2,sp2,vx2,vp2,mchi,lam,nu,gg)

    sp2 = -sp2
    vp2 = -vp2

    print('@',format(istep,'10'),format(sp0[3,3]-sp2[3,3],' .16E'),\
            format(sx0[15,2] - sx2[15,2],' .16E'),\
            format(vp0[2,10,0] - vp2[2,10,0],' .16E'),\
            format(vp0[8,4,1] - vp2[8,4,1],' .16E'),\
            format(vx0[4,5,0] - vx2[4,5,0],' .16E'),\
            format(vx0[12,7,1] - vx2[12,7,1],' .16E'))
#

    h1 = hamil(sx1,sp1,vx1,vp1,mchi,lam,nu,gg)
    # print('{:.16e}'.format(h1),'{:.16e}'.format(h0),'{:.16e}'.format(h1-h0))

    rho = min(1.0,np.exp(h0-h1))
    rand_num = np.random.random()

    itry = itry + 1
    if (rand_num <= rho):
        iacc = iacc + 1
    else:
        sx1 = sx0
        vx1 = vx1

    if (istep >= NTHERM and mod(istep,NSKIP) == 0):
        print(format(istep,'10'),format(sx1[8,8],' :16E'))

    sx0 = sx1
    vx0 = vx1

pacc = float(iacc/itry)*100
print('# HMC Metropolis test statistics')
print('# tau = %5d, nmd = %5d, dt = %.16E'%(tau,nmd,tau/nmd))
print('# iacc = %10d, ityr = %10d, Pacc = %5d'%(iacc,itry,pacc))
print('# mchi = %.16E, lam = %.16E'%(mchi,lam))
print('#  nu  = %.16E,  g  = %.16E'%(nu,gg))
