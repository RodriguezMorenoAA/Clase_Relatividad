# %%
%pip install -q --upgrade pymcel rebound montu

# %% [markdown]
# # Aberración de la luz en estrellas

# %%
import numpy as np
import matplotlib.pyplot as plt
import montu as mn
import rebound as rb
import pymcel as pc

# %%
tabla, jd, X = pc.consulta_horizons(id='399',
                                    location='@0',
                                    epochs ='2025-03-24 15:40:00'
                )


# %%
v_tierra = X[3:]
v_tierra

# %% [markdown]
# Ahora el beta de la estrella

# %%
beta_vec = -v_tierra / pc.constantes.c
beta = np.linalg.norm(beta_vec)
beta_vec, beta 

# %%
gamma = 1/np.sqrt(1-beta**2)
gamma

# %%
allstars = mn.Stars()

# %%
star = allstars.get_stars(ProperName='Aldebaran')
star

# %%
ra = np.array(star.data.RAJ2000)[0]
dec = np.array(star.data.DecJ2000)[0]
ra, dec

# %%
import spiceypy as spy
deg = np.pi/180
rad = 1 / deg

# %%
nprima_equ = spy.latrec(1, ra*15*deg, dec*deg)
nprima_equ

# %% [markdown]
# ### Transformar al sistema eclíptico de coordenadas

# %%
Requ2ecl = spy.pxform('J2000', 'ECLIPJ2000', 0) #Matriz de rotación
Requ2ecl

# %% [markdown]
# Ya se puede convertir el vector

# %%
nprima_ecl = spy.mxv(Requ2ecl, nprima_equ)
nprima_ecl
n = (nprima_equ + (((gamma -1)/beta**2)*np.dot(beta_vec, nprima_equ) + gamma)*beta_vec) / (gamma*(1+np.dot(beta_vec, nprima_equ)))

# %% [markdown]
# ### Quiero calcular la ascención recta y la declinación después de la aberración

# %%
n_equ = spy.mxv(spy.invert(Requ2ecl),n)
n_equ

# %%
r, long, lat = spy.reclat(n_equ)

# %%
ra_aberrada = long * rad / 15
dec_aberrada = lat * rad
mn.Util.dec2hex(ra_aberrada), mn.Util.dec2hex(dec_aberrada)

# %%


# %%



