import numpy as np
import pb
m, P, cnt = pb.source_line_means(pb.H1_SOURCES + ['H1'], pb.panel_genes())
np.savez(pb.CACHE / 'line_means_abcsources.npz', mean=m, per=P, cnt=cnt, sources=np.array(pb.H1_SOURCES + ['H1']))
print('abc extension genes (nonzero):', int((m != 0).sum()))
