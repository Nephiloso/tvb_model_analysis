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
idx = np.where((mapping>72)&(mapping<75))[0]
new_coupling = {}
new_coupling['__header__'] = connect['__header__']
new_coupling['__version__'] = connect['__version__']
new_coupling['__globals__'] = connect['__globals__']
new_coupling['LocalCoupling'] = coupling[idx,:][:,idx]
savemat("local_connectivity_843.mat", new_coupling)
new_connect = loadmat('./local_connectivity_843.mat')

triang = np.loadtxt(triang_path,dtype=np.int32) # 不用int32下面生成vertex_triangles会出错
vert = np.loadtxt(vert_path)  # * 先生成vertex_triangles，不包含上述idx的就排除
vert_n = np.loadtxt(vert_n_path)  # * 先生成vertex_triangles，不包含上述idx的就排除
vertex_triangles0 = [[] for _ in range(vert.shape[0])]
for k in range(triang.shape[0]):  #function vertex_triangles, _find_vertex_triangles in surface: map the vertex onto the surface - assign the index of the triangel to the vertex
    vertex_triangles0[triang[k, 0]].append(k)
    vertex_triangles0[triang[k, 1]].append(k)
    vertex_triangles0[triang[k, 2]].append(k)
vertex_triangles = np.array([np.array(xi) for xi in vertex_triangles0])
# * 根据idx生成新的map和vert_tri
new_map = mapping[idx] - 73
new_vertex_triangles = vertex_triangles[idx] # 为什么有的元素长度>3?说明有很多triangle共享这个vertex
new_vert = vert[idx]
new_vert_n = vert_n[idx]
print(new_vertex_triangles.shape)
tri_idx = np.unique(np.concatenate(new_vertex_triangles))
# * 找到idx对应的vertex_triangle，保留相应位置的triang
new_triang = triang[tri_idx]

# save the new data
np.savetxt('regionMapping_843_2.txt', new_map, fmt='%.0f', delimiter=' ')
np.savetxt('triangles0.txt', new_triang,fmt='%d')
np.savetxt('vertices0.txt', new_vert, fmt='%.6f')
np.savetxt('vertex_normals0.txt', new_vert_n, fmt='%.6f')

