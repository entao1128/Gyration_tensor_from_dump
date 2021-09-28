#!/usr/bin/env python
# coding: utf-8

# In[515]:


import numpy as np
import MD_IO_tool as tool
import numpy.linalg as la



def extract_polymer(atom):
    # atom save the particle information of the configuration, where atom[i, :] is the corresponding info for 
    # particle i, atom[:, 0] is mol_id, atom[:, 1] is type, atom[:, 2:] are position
    # polymers are defined as type_id > 1
    
    # extract the polymer chains, num of chains, chain length
    polymer = atom[atom[:, 0]>1]
    chain_num = int(max(atom[:,0]) - 1)
    chl_one = atom[atom[:, 0] == 2]   # chain length
    chl = np.size(chl_one[:, 0])
    polymer = polymer.reshape(chain_num, chl, 5)
    return polymer, chain_num, chl  # return polymers, num of chains, chain length


# In[334]:


def orientation_analysis(atom, cutoff, radius, boxh):
    # atom save the particle information of the configuration, where atom[i, :] is the corresponding info for 
    # particle i, atom[:, 0] is mol_id, atom[:, 1] is type, atom[:, 2:] are position
    # cutoff is the distance for bound layer, radius is the NP raidus
    
    # extract polymer infomation, return polymer[chain_id, chain_order, info]
    polymer, chain_num, chl = extract_polymer(atom)
     
    # extract bound and bulk polymer
    bound_poly, bulk_poly = extract_bound(polymer, chain_num, chl, cutoff, radius, boxh)
    #bound_poly, bulk_poly = extract_bound2(polymer, chain_num, chl, cutoff, radius)

    # measure number of bulk polymer chains and bound polymer chains
    n_bulk = np.size(bulk_poly[:, 0, 0])
    n_bound = np.size(bound_poly[:, 0, 0])
    
    # calculate the center of mass for bound and bulk polymer 
    bound_center = calc_center(bound_poly, boxh)
    bulk_center = calc_center(bulk_poly, boxh)
    
    #bound_poly, bulk_poly = extract_bound(polymer, chain_num, chl, cutoff, radius, boxh)
    print(bound_poly.shape)
    print(bulk_poly.shape)
    
    # substract the center of mass for bound and bulk polymer
    for i in range(n_bound):
        for j in range(3):
            for k in range(chl):
                dr = bound_poly[i, k, j+2] - bound_center[i, j]
                if (dr > boxh):
                    bound_poly[i, k, j+2] -= (boxh * 2)
                elif (dr < -boxh):
                    bound_poly[i, k, j+2] += (boxh * 2)
                bound_poly[i, k, j+2] -= bound_center[i, j]
            
    for i in range(n_bulk):
        for j in range(3):
            for k in range(chl):
                dr = bulk_poly[i, k, j+2] - bulk_center[i, j]
                if (dr > boxh):
                    bulk_poly[i, k, j+2] -= (boxh * 2)
                elif (dr < -boxh):
                    bulk_poly[i, k, j+2] += (boxh * 2)
                bulk_poly[i, k, j+2] -= bulk_center[i, j]

    # calculate the gyration tensor for bound and bulk polymer
    bound_gy_tensor = calc_gyration_tensor(bound_poly)
    bulk_gy_tensor = calc_gyration_tensor(bulk_poly)
    
    # calculate eignevalue and eigenvectors
    bound_eigvecs, bound_eigen_id = calc_eigen(bound_gy_tensor)
    bulk_eigvecs, bulk_eigen_id = calc_eigen(bulk_gy_tensor)
    
    # Calculate gyration tensor parameters: Rg2, eta, c, k2
    # eta: asphericity, c: acylindricity, k2: anisotropy
    bound_Rg2, bound_para = calc_gytensor_para(bound_gy_tensor)
    bulk_Rg2, bulk_para = calc_gytensor_para(bulk_gy_tensor)
    
    # calculate the eigenvalues of gyration tensor
    bound_eigvals = np.zeros((n_bound, 3))
    bulk_eigvals = np.zeros((n_bulk, 3))
    for i in range(n_bound):
        bound_eigvals[i, :] = list(sorted(la.eigvalsh(bound_gy_tensor[:, :, i]), reverse=True))
    for i in range(n_bulk):
        bulk_eigvals[i, :] = list(sorted(la.eigvalsh(bulk_gy_tensor[:, :, i]), reverse=True))
    #print(bound_eigvals)

    # calculate the orientation order parameter
    bound_orient = calc_orient_order(bound_center, bound_eigvecs, bound_eigen_id)
    bulk_orient = calc_orient_order(bulk_center, bulk_eigvecs, bulk_eigen_id)

    #Ly = np.sqrt()
    #cos_angle = bound_center[i, 0
    #bound_theta[i] = 
    return n_bound, n_bulk, bound_Rg2, bulk_Rg2, bound_para, bulk_para, bound_orient, bulk_orient, bound_eigvals, bulk_eigvals
    #return n_bound, n_bulk, bound_center, bulk_center, bound_para, bulk_para, bound_orient, bulk_orient


# In[251]:


def calc_eigen(tensor):
    n_poly = np.size(tensor[0,0, :])
    eigvecs = np.zeros((n_poly, 4, 3))
    eigen_id = np.zeros((n_poly, 3))
    for i in range(n_poly):
        eigvecs[i, 0, :], eigvecs[i, 1:, :] = la.eig(tensor[:, :, i])
        #bound_eigvecs[i, 0, :] = list(sorted(bound_eigvecs[:, 0, :], reverse=True))
        eigen_id[i, :] = np.argsort(-eigvecs[i, 0, :])  # sort eigenvalues descending
    return eigvecs, eigen_id


# In[260]:


def calc_gyration_tensor(poly):
    # calculate the gyration tensor for the target polymer group
    
    # check numer of chains and chain length
    n_poly = np.size(poly[:, 0, 0])
    chl = np.size(poly[0, :, 0])
    
    # define the gyration tensor
    gy_tensor = np.zeros((3, 3, n_poly))
    
    #calculate gyration tensor
    count = 0
    for n in range(n_poly):
        item = poly[n, :, :]
        for i in range(3):
            for j in range(i, 3):
                gy_tensor[i, j, count] = np.dot(item[:, i+2], item[:, j+2])
                gy_tensor[j, i, count] = gy_tensor[i, j, count]
        count +=1 
    gy_tensor /= chl  
    
    return gy_tensor


# In[310]:


def extract_bound2(polymer, chain_num, chl, cutoff, radius):
    # define bound layer based on distance to NP surface
    
    # initialize list to save bound layer and bulk polymers
    bound_poly = []
    bulk_poly = []
    
    for i in range(0, chain_num):
        for j in range(0, chl):
            dr2 = polymer[i, j, 2]*polymer[i, j, 2]+polymer[i, j, 3]*polymer[i, j, 3]+polymer[i, j, 4]*polymer[i, j, 4]
            dr = np.sqrt(dr2)
            if (dr <= cutoff+radius):
                bound_poly.append(polymer[i, :, :])
                break
            elif (j == chl-1):
                bulk_poly.append(polymer[i, :, :])
    
    # change the list to array for post analysis
    bulk_poly = np.array(bulk_poly)
    bound_poly = np.array(bound_poly)
    
    return bound_poly, bulk_poly


# In[335]:


def extract_bound(polymer, chain_num, chl, cutoff, radius, boxh):
    # define bound layer based on distance to NP surface
    
    # initialize list to save bound layer and bulk polymers
    bound_poly = []
    bulk_poly = []
    # calculate center of mass

    center = calc_center(polymer, boxh)

    for i in range(0, chain_num):
        for j in range(0, chl):
            dr2 = polymer[i, j, 2]*polymer[i, j, 2]+polymer[i, j, 3]*polymer[i, j, 3]+polymer[i, j, 4]*polymer[i, j, 4]
            dr = np.sqrt(dr2)           
            if (dr <= cutoff+radius):
                bound_poly.append(polymer[i, :, :])
                break
            elif (j == chl-1):
                dr2 = center[i, 0]*center[i,0] + center[i,1]*center[i,1] + center[i,2]*center[i,2]
                dr = np.sqrt(dr2)
                if (dr <= cutoff+radius+8.0):
                    bulk_poly.append(polymer[i, :, :])
    
    # change the list to array for post analysis
    bulk_poly = np.array(bulk_poly)
    bound_poly = np.array(bound_poly)
    
    return bound_poly, bulk_poly


# In[304]:


def calc_center(poly, boxh):
    # calculate the center of mass of each chain
    # poly save polymer information of a certain group    
    n_poly = np.size(poly[:, 0, 0])
    chl = np.size(poly[0, :, 0])
    
    polymer = np.zeros((n_poly, chl, 3))
    # calculate the center of mass for the target group
    center = np.zeros((n_poly, 3))    
    for i in range(n_poly):
        for j in range(3):
            for k in range(chl):
                if (k==0):
                    center[i, j] = poly[i, 0, j+2]
                else:
                    dr = poly[i, k, j+2] - center[i, j]
                    if (dr > boxh):
                        polymer[i, k, j] = poly[i, k, j+2] - boxh*2
                    elif (dr < -boxh):
                        polymer[i, k, j] = poly[i, k, j+2] + boxh*2
                    else:
                        polymer[i, k, j] = poly[i, k, j+2]
            center[i, j] = np.mean(polymer[i, :, j])
            if (center[i, j] > boxh):
                center[i, j] -= (2*boxh)
            elif (center[i, j] < -boxh):
                center[i, j] += (2*boxh)
    return center


# In[9]:


def calc_gytensor_para(tensor):
    # Calculate gyration tensor parameters: Rg2, eta, c, k
    
    # check the chain number
    n_poly = np.size(tensor[0,0, :])

    # calculate the eigenvalues of gyration tensor (this way is more efficient, for eigenvalues only)
    # ordered descending
    eigvals = np.zeros((n_poly, 3))
    for i in range(n_poly):
        eigvals[i, :] = list(sorted(la.eigvalsh(tensor[:, :, i]), reverse=True))
    
    # Rg2: squared gyration radius 
    Rg2 = np.zeros(n_poly)
    for i in range(n_poly):
        Rg2[i] = sum(eigvals[i,:])
        
    # eta: asphericity, c: acylindricity, k2: anisotropy
    tensor_para = np.zeros((n_poly, 3))
    for i in range(n_poly):
        if (Rg2[i]>0):
            tensor_para[i, 0] = (eigvals[i,0] - 0.5 * (eigvals[i,1]+eigvals[i,2])) / Rg2[i] # eta
            tensor_para[i, 1] = (eigvals[i,1]-eigvals[i,2]) / Rg2[i]                # c
            tensor_para[i, 2] = (tensor_para[i, 0]**2 + 0.75*tensor_para[i, 1]**2)                # k2
    
    return Rg2, tensor_para


# In[10]:


def calc_orient_order(center, eigvecs, eigen_id):
    # calculate the orientation order parameter
    
    # check num of polymer chains
    n_poly = np.size(center[:, 0])
    
    theta = np.zeros((n_poly, 3))
    for i in range(n_poly):
        Lx = np.sqrt(center[i,:].dot(center[i,:]))
        for j in range(3):
            e_id = int(eigen_id[i, j])
            ei = eigvecs[i, e_id+1, :]
            Ly = np.sqrt(ei.dot(ei))
            cos_angle = center[i,:].dot(ei)/(Lx*Ly)
            theta[i, e_id] = np.arccos(cos_angle)

    orient_para = np.zeros((n_poly, 3))
    for i in range(n_poly):
        for j in range(3):
            orient_para[i, j] = 0.5 * (3 * (np.cos(theta[i, j])*np.cos(theta[i, j]))-1)
    
    mean_orient = np.zeros(3)
    for i in range(3):
        mean_orient[i] = np.mean(orient_para[:, i])  
    
    return mean_orient


