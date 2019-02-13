import numpy as np
########### constant
pi = np.pi
NSITE = 16
mchi = np.sqrt(5)
gg = 0.8
lam = 0.1
nu = 3.0
NTHERM = 100
NSKIP = 1
NSAMPLE = 0
nmd = 10
seed = 98

########### set up of unit vector and periodic vector
mu = np.array([[1,0],[0,1]],dtype='int64')

pn = [0] * NSITE
for i in range(NSITE): # i = 0...15
    pn[i] = i
pn = np.insert(pn,NSITE,[0,1])
pn = np.insert(pn,0,[NSITE-2,NSITE-1])

################################################
# computation of hamiltonian
################################################
def hamil (sx,sp,vx,vp,mchi,lam,nu,gg) :
    h = 0.0
    for j in range(NSITE) :
        for i in range(NSITE) :
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

    update_x(dt2,sx,sp,vx,vp,mchi,lam,nu,gg)

    for i in range(nmd-1):
        update_p(dt,sx,sp,vx,vp,mchi,lam,nu,gg)
        update_x(dt,sx,sp,vx,vp,mchi,lam,nu,gg)

    update_p(dt,sx,sp,vx,vp,mchi,lam,nu,gg)
    update_x(dt2,sx,sp,vx,vp,mchi,lam,nu,gg)


def update_x (dt,sx,sp,vx,vp,mchi,lam,nu,gg):
    sx = sx + sp*dt
    vx = vx + vp*dt


def update_p (dt,sx,sp,vx,vp,mchi,lam,nu,gg):
    sf = get_scalar_force(sx,vx,mchi,lam,nu,gg)
    sp = sp + sf*dt

    vf = get_vector_force(sx,vx,mchi,lam,nu,gg)
    vp = vp + vf*dt


######################################################
###### calculation of force
######################################################
def get_scalar_force (sx,vx,mchi,lam,nu,gg):
    sf = [0]*NSITE*NSITE
    sf = np.reshape(sf,(NSITE,NSITE),order='C')

    for j in range(NSITE):
        for i in range(NSITE):
            div_s = div_sca(sx,i,j)
            div_v = div_vec(vx,i,j)
            sf[i,j] = div_s - 4.0*lam*sx[i,j]*(sx[i,j]**2 - nu**2) + gg*(div_v - gg*sx[i,j])/mchi**2
    return sf


def get_vector_force (sx,vx,mchi,lam,nu,gg):
    vf = [0]*NSITE*NSITE*2
    vf = np.reshape(vf,(NSITE,NSITE,2),order='C')

    for k in range(2):
        for j in range(NSITE):
            for i in range(NSITE):
                pni = pn[i + mu[k,0]]
                pnj = pn[j + mu[k,1]]
                div_v = div_vec(vx,pni,pnj)

                vf[i,j,k] = -vx[i,j,k] + 0.5*(div_v - gg*sx[pni,pnj])/mchi**2

                pni = pn[i - mu[k,0]]
                pnj = pn[j - mu[k,1]]
                div_v = div_vec(vx,pni,pnj)

                vf[i,j,k] = +vf[i,j,k] - 0.5*(div_v - gg*sx[pni,pnj])/mchi**2

    return vf


######################################################
####### derivative part
######################################################
def div_sca (s,i,j) :
    div_s = 0.0
    for k in range(2) :
        div_s = div_s + s[pn[i+mu[k,0]],pn[j+mu[k,1]]] - 2.0*s[i,j] + s[pn[i-mu[k,0]],pn[j-mu[k,1]]]
    return div_s


def div_vec (v,i,j) :
    div_v = 0.0
    for k in range(2) :
        div_v = div_v + v[pn[i+mu[k,0]],pn[j+mu[k,1]],k] - v[pn[i-mu[k,0]],pn[j-mu[k,1]],k]

    return div_v


def sum_vec (v,i,j) :
    sum_v = 0.0
    for k in range(2) :
        sum_v = sum_v + v[i,j,k]**2
    return sum_v


################################################
################################################
################ main program ##################
################################################
################################################

sx0 = [0]*NSITE*NSITE
sx0 = np.random.normal(0,1,NSITE*NSITE)
sx0 = np.reshape(sx0,(NSITE,NSITE),order='C')

vx0 = [0]*NSITE*NSITE*2
vx0 = np.random.normal(0,1,NSITE*NSITE*2)
vx0 = np.reshape(vx0,(NSITE,NSITE,2),order='C')

np.random.seed(seed)
tau = 1.0
itry = 0
iacc = 0

for istep in range(NTHERM + NSKIP*NSAMPLE):
    sp0 = [0]*NSITE*NSITE
    sp0 = np.random.normal(0,1,NSITE*NSITE)
    sp0 = np.reshape(sp0,(NSITE,NSITE),order='C')

    vp0 = [0]*NSITE*NSITE*2
    vp0 = np.random.normal(0,1,NSITE*NSITE*2)
    vp0 = np.reshape(vp0,(NSITE,NSITE,2),order='C')

    h0 = hamil(sx0,sp0,vx0,vp0,mchi,lam,nu,gg)
    # print(h0)

    sx1 = sx0
    sp1 = sp0
    vx1 = vx0
    vp1 = vp0

    leapfrog_xpx_md(tau,nmd,sx1,sp1,vx1,vp1,mchi,lam,nu,gg)

    sp1 = -sp1
    vp1 = -vp1

# _CHECK_REVERSE_
    # sx2 = sx1
    # sp2 = sp1
    # vx2 = vx1
    # vp2 = vp1
    #
    # leapfrog_xpx_md(tau,nmd,sx2,sp2,vx2,vp2,mchi,lam,nu,gg)
    # h2 = hamil(sx2,sp2,vx2,vp2,mchi,lam,nu,gg)
    #
    # sp2 = -sp2
    # vp2 = -vp2
    #
    # print('@',istep,sp0[3,3]-sp2[3,3],sx0[15,2]-sx2[15,2],\
    #         vp0[2,10,0] - vp2[2,10,0], vp0[8,4,1] - vp2[8,4,1],\
    #         vx0[4,5,0] - vx2[4,5,0], vx0[12,7,1] - vx1[12,7,1])
#

    h1 = hamil(sx1,sp1,vx1,vp1,mchi,lam,nu,gg)
    print('{:.16e}'.format(h1),'{:.16e}'.format(h0),'{:.16e}'.format(h1-h0))
    rho = min(1.0,np.exp(h0-h1))
    rand_num = np.random.random()

    itry = itry + 1
    if (rand_num <= rho):
        iacc = iacc + 1
    else:
        sx1 = sx0
        vx1 = vx1

    sx0 = sx1
    vx0 = vx1

pacc = float(iacc/itry)*100
print(iacc,itry,'{:.15f}'.format(pacc))
