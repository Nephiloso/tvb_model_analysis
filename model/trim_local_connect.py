import numpy as np
import os
from scipy.io import loadmat,savemat

os.chdir(r'C:\Users\wpp_1\Documents\Neurasmus\VU\Internship\codes')
mapping_path = './regionMapping_16k_76.txt'
triang_path = './triangles.txt'
vert_path = './vertices.txt'
vert_n_path = './vertex_normals.txt'
connect = loadmat('./local_connectivity_16384.mat')
coupling = connect['LocalCoupling']
mapping = np.loadtxt(mapping_path) # * 找到73 74的idx
idx = np.where((mapping>72)&(mapping<74))[0]
new_coupling = {}
new_coupling['__header__'] = connect['__header__']
new_coupling['__version__'] = connect['__version__']
new_coupling['__globals__'] = connect['__globals__']
new_coupling['LocalCoupling'] = coupling[idx,:][:,idx]
savemat("local_connectivity_180.mat", new_coupling)
new_connect = loadmat('./local_connectivity_180.mat')

new_map = mapping[idx] - 73

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
new_vert = vert[idx]
new_vert_n = vert_n[idx]
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
new_triang = np.delete(new_triang0,zero_idx, axis=0)

# save the new data
np.savetxt('regionMapping_180_1.txt', new_map, fmt='%.0f', delimiter=' ')
np.savetxt('triangles0.txt', new_triang,fmt='%d')
np.savetxt('vertices0.txt', new_vert, fmt='%.6f')
np.savetxt('normals0.txt', new_vert_n, fmt='%.6f')  # two kinds of normals of Surface in TVB: vertex normals -> normals; trialgel normals are computed by TVB