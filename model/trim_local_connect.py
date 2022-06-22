import numpy as np
import os
from scipy.io import loadmat,savemat
from scipy.sparse import csr_matrix

def delete_from_csr(mat, row_indices=[], col_indices=[]):
    """
    Remove the rows (denoted by ``row_indices``) and columns (denoted by ``col_indices``) from the CSR sparse matrix ``mat``.
    """
    if not isinstance(mat, csr_matrix):
        raise ValueError("works only for CSR format -- use .tocsr() first")

    rows = []
    cols = []
    if row_indices:
        rows = list(row_indices)
    if col_indices:
        cols = list(col_indices)

    if len(rows) > 0 and len(cols) > 0:
        row_mask = np.ones(mat.shape[0], dtype=bool)
        row_mask[rows] = False
        col_mask = np.ones(mat.shape[1], dtype=bool)
        col_mask[cols] = False
        return mat[row_mask][:,col_mask]
    elif len(rows) > 0:
        mask = np.ones(mat.shape[0], dtype=bool)
        mask[rows] = False
        return mat[mask]
    elif len(cols) > 0:
        mask = np.ones(mat.shape[1], dtype=bool)
        mask[cols] = False
        return mat[:,mask]
    else:
        return mat

os.chdir(r'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes')
mapping_path = './regionMapping_16k_76.txt'
triang_path = './triangles.txt'
vert_path = './vertices.txt'
vert_n_path = './vertex_normals.txt'
connect = loadmat('./local_connectivity_16384.mat')
coupling = connect['LocalCoupling']
mapping = np.loadtxt(mapping_path)
idx = np.where((mapping>72)&(mapping<74))[0] # * 找到73 74的idx

new_map0 = mapping[idx] - 73

triang = np.loadtxt(triang_path,dtype=np.int32) # 不用int32下面生成vertex_triangles会出错
vert = np.loadtxt(vert_path)  # * 先生成vertex_triangles，不包含上述idx的就排除
vert_n = np.loadtxt(vert_n_path)  # * 先生成vertex_triangles，不包含上述idx的就排除
vertex_triangles0 = [[] for _ in range(vert.shape[0])]

# Get the verteices that are used in the triangles
for k in range(triang.shape[0]):  #function vertex_triangles, _find_vertex_triangles in surface: map the vertex onto the surface - assign the index of the triangel to the vertex
    vertex_triangles0[triang[k, 0]].append(k)
    vertex_triangles0[triang[k, 1]].append(k)
    vertex_triangles0[triang[k, 2]].append(k)
vertex_triangles = np.array([np.array(xi) for xi in vertex_triangles0],dtype=object)
# * 根据idx生成新的map和vert_tri
new_vertex_triangles = vertex_triangles[idx] # 为什么有的元素长度>3?说明有很多triangle共享这个vertex
new_vert0 = vert[idx]
new_vert_n0 = vert_n[idx]
tri_idx = np.unique(np.concatenate(new_vertex_triangles))
# * 找到idx对应的vertex_triangle，保留相应位置的triang
mid_triang = triang[tri_idx]
new_triang0 = np.zeros(mid_triang.shape)

# a. weive the triangels containing vertex that is out of the region b. remap the idx of triangels
for l in range(mid_triang.shape[0]):
    try:
        new_triang0[l,0] = np.argwhere(idx==mid_triang[l][0])[0][0]
        new_triang0[l,1] = np.argwhere(idx==mid_triang[l][1])[0][0]
        new_triang0[l,2] = np.argwhere(idx==mid_triang[l][2])[0][0]
    except IndexError:
        new_triang0[l,:] = [0,0,0]
zero_idx = np.argwhere(np.all(new_triang0[:,...]==0,axis=1))
new_triang1 = np.delete(new_triang0,zero_idx, axis=0)

# Delete isolated triangles
iso_tris = []
for i in range(len(idx)):
    posit = np.where(new_triang1==i)[0]
    if len(posit) == 1:
        iso_tris.append(posit[0])
iso_tris_u = np.unique(iso_tris)
isos = []
for iso_tri in iso_tris_u:
    if np.count_nonzero(iso_tris==iso_tri)==3:
        isos.append(new_triang1[iso_tri][0])
        isos.append(new_triang1[iso_tri][1])
        isos.append(new_triang1[iso_tri][2])
        new_triang1 = np.delete(new_triang1, iso_tri,axis=0)
isos=np.array(isos,dtype=int)
new_vert = np.delete(new_vert0,isos,axis=0)
new_vert_n = np.delete(new_vert_n0,isos,axis=0)
new_map = np.delete(new_map0,isos,axis=0)

idx_mapping1 = np.array(range(new_vert0.shape[0]))
idx_mapping1 = np.delete(idx_mapping1,isos,axis=0)
new_triang = np.zeros_like(new_triang1)
# remap the idx for the triangle again
for l in range(new_triang1.shape[0]):
    new_triang[l,0] = np.argwhere(idx_mapping1==new_triang1[l][0])[0][0]
    new_triang[l,1] = np.argwhere(idx_mapping1==new_triang1[l][1])[0][0]
    new_triang[l,2] = np.argwhere(idx_mapping1==new_triang1[l][2])[0][0]

# Get the local coupling matrix
nc0 = coupling[idx,:][:,idx]
nc = delete_from_csr(nc0.tocsr(),isos.tolist(),isos.tolist()).tocsr()

new_coupling = {}
new_coupling['__header__'] = connect['__header__']
new_coupling['__version__'] = connect['__version__']
new_coupling['__globals__'] = connect['__globals__']
new_coupling['LocalCoupling'] = nc

# save the new data
savemat("local_connectivity_177.mat", new_coupling)
np.savetxt('regionMapping_177_1.txt', new_map, fmt='%.0f', delimiter=' ')
np.savetxt('triangles0.txt', new_triang,fmt='%d')
np.savetxt('vertices0.txt', new_vert, fmt='%.6f')
np.savetxt('normals0.txt', new_vert_n, fmt='%.6f')  # two kinds of normals of Surface in TVB: vertex normals -> normals; trialgel normals are computed by TVB