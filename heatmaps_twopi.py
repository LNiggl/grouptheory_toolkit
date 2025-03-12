import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pylab as plt

if __name__ == '__main__':
    generic_label_string  =     [
                                r'$\hat{a}_1^{\,u}$',
                                r'$\hat{a}_1^{\,g}$',
                                r'$\hat{x}^{\,u}$',r'$\hat{y}^{\,u}$',r'$\hat{z}^{\,u}$',
                                r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                r'$\hat{a}_2^{\,u}$',
                                r'$\hat{a}_2^{\,g}$',
                                r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
                                r'$\hat{\epsilon}_1^{\,u}$',r'$\hat{\epsilon}_2^{\,u}$',
                                r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                ]
    ##  1 0 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi100_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi100_vecdata_irreps_seq.npy")
    # print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")

    labels_v = [r'$\hat{a}_1^{\,g}$',r'$\hat{x}^{\,u}$',r'$\hat{y}^{\,u}$',r'$\hat{z}^{\,u}$',r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(1,0,0)}' for x in range(len(vecs))]
    labels_irreps = [r'$A_1^g$',r'$T_1^u$',r'$E^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y],rotation = 0)

    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5, 5], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)

    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,5.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi100_heatmap.png")
    
    ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(1,0,0)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi100_heatmap_yaxis_idx.png")

    plt.close()

    #  1 1 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi110_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi110_vecdata_irreps_seq.npy")
    print("irreps 110:",irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u}$',r'$\hat{y}^{\,u}$',r'$\hat{z}^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(1,1,0)}' for x in range(len(vecs))]
    labels_irreps = [r'$A_1^g$',r'$T_1^u$',r'$T_2^u$', r'$T_2^g$',r'$E^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y], rotation = 0)
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,11], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi110_heatmap.png")
    
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(1,1,0)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi110_heatmap_yaxis_idx.png")
    plt.close()

    #  1 1 1  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi111_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi111_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u}$',r'$\hat{y}^{\,u}$',r'$\hat{z}^{\,u}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
    ]
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(1,1,1)}' for x in range(len(vecs))]
    labels_irreps = [r'$A_1^g$',r'$T_1^u$',r'$A_2^u$', r'$T_2^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y], rotation = 0)
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5, 4.5, 6.5], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,5,7.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi111_heatmap.png")
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(1,1,1)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi111_heatmap_yaxis_idx.png")
    
    plt.close()

    #  2 1 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi210_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi210_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,g}$',
                                        r'$\hat{\tau}_1^{\,u,(1)}$',r'$\hat{\tau}_2^{\,u,(1)}$',r'$\hat{\tau}_3^{\,u,(1)}$',
                                        r'$\hat{\tau}_1^{\,u,(2)}$',r'$\hat{\tau}_2^{\,u,(2)}$',r'$\hat{\tau}_3^{\,u,(2)}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
                                        r'$\hat{\epsilon}_1^{\,g,(1)}$',r'$\hat{\epsilon}_2^{\,g,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(2)}$',r'$\hat{\epsilon}_2^{\,g,(2)}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^g$', r'$T_2^{u,(1)}$', r'$T_2^{u,(2)}$', r'$T_2^g$', r'$E^{g,(1)}$',  r'$E^{g,(2)}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(2,1,0)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.075)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 53, width = 1)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi210_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    # ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_yticklabels(["1","","","","","6","","","","","","12","","","","","","18",
                        "","","","","","24"], size = 15)
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(2,1,0)}$", size = 15)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi210_heatmap_yaxis_idx.png")
    
    plt.close()

    ##  2 1 1  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi211_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi211_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$',                                        
                                        r'$\hat{\epsilon}_1^{\,u}$',r'$\hat{\epsilon}_2^{\,u}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^u$', r'$T_2^{u}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$', r'$E^{u}$',  r'$E^{g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(2,1,1)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.075)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 53, width = 1)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi211_heatmap.png",dpi = 100)
    
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(2,1,1)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi211_heatmap_yaxis_idx.png")
    plt.close()

    # 2 2 1 ##

    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi221_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi221_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$',                                        
                                        r'$\hat{\epsilon}_1^{\,u}$',r'$\hat{\epsilon}_2^{\,u}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^u$', r'$T_2^{u}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$', r'$E^{u}$',  r'$E^{g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(2,2,1)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.075)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 53, width = 1)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi221_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(2,2,1)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi221_heatmap_yaxis_idx.png")
    
    plt.close()
    
    # 3 1 0 ## 

    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi310_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi310_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,g}$',
                                        r'$\hat{\tau}_1^{\,u,(1)}$',r'$\hat{\tau}_2^{\,u,(1)}$',r'$\hat{\tau}_3^{\,u,(1)}$',
                                        r'$\hat{\tau}_1^{\,u,(2)}$',r'$\hat{\tau}_2^{\,u,(2)}$',r'$\hat{\tau}_3^{\,u,(2)}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
                                        r'$\hat{\epsilon}_1^{\,g,(1)}$',r'$\hat{\epsilon}_2^{\,g,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(2)}$',r'$\hat{\epsilon}_2^{\,g,(2)}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^g$', r'$T_2^{u,(1)}$', r'$T_2^{u,(2)}$', r'$T_2^g$', r'$E^{g,(1)}$',  r'$E^{g,(2)}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(3,1,0)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.075)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 53, width = 1)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi310_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(3,1,0)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi310_heatmap_yaxis_idx.png")
    
    plt.close()

    ##  3 1 1  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi311_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi311_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$',                                        
                                        r'$\hat{\epsilon}_1^{\,u}$',r'$\hat{\epsilon}_2^{\,u}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^u$', r'$T_2^{u}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$', r'$E^{u}$',  r'$E^{g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(3,1,1)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.075)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 53, width = 1)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi311_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(3,1,1)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi311_heatmap_yaxis_idx.png")
    
    plt.close()

    # 3 2 1 ## 
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi321_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi321_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,u}$',
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,u,(3)}$',r'$\hat{y}^{\,u,(3)}$',r'$\hat{z}^{\,u,(3)}$',
                                        r'$\hat{x}^{\,g,(1)}$',r'$\hat{y}^{\,g,(1)}$',r'$\hat{z}^{\,g,(1)}$',
                                        r'$\hat{x}^{\,g,(2)}$',r'$\hat{y}^{\,g,(2)}$',r'$\hat{z}^{\,g,(2)}$',
                                        r'$\hat{x}^{\,g,(3)}$',r'$\hat{y}^{\,g,(3)}$',r'$\hat{z}^{\,g,(3)}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{a}_2^{\,g}$',
                                        r'$\hat{\tau}_1^{\,u,(1)}$',r'$\hat{\tau}_2^{\,u,(1)}$',r'$\hat{\tau}_3^{\,u,(1)}$',
                                        r'$\hat{\tau}_1^{\,u,(2)}$',r'$\hat{\tau}_2^{\,u,(2)}$',r'$\hat{\tau}_3^{\,u,(2)}$',
                                        r'$\hat{\tau}_1^{\,u,(3)}$',r'$\hat{\tau}_2^{\,u,(3)}$',r'$\hat{\tau}_3^{\,u,(3)}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$', 
                                        r'$\hat{\tau}_1^{\,g,(3)}$',r'$\hat{\tau}_2^{\,g,(3)}$',r'$\hat{\tau}_3^{\,g,(3)}$',                                       
                                        r'$\hat{\epsilon}_1^{\,u,(1)}$',r'$\hat{\epsilon}_2^{\,u,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,u,(2)}$',r'$\hat{\epsilon}_2^{\,u,(2)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(1)}$',r'$\hat{\epsilon}_2^{\,g,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(2)}$',r'$\hat{\epsilon}_2^{\,g,(2)}$'
                                        ]
    labels_irreps = [r'$A_1^u$',r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,  r'$T_1^{u,(3)}$', r'$T_1^{g,(1)}$', r'$T_1^{g,(2)}$' ,  r'$T_1^{g,(3)}$',
                     r'$A_2^u$', r'$A_2^g$', r'$T_2^{u,(1)}$', r'$T_2^{u,(2)}$' ,  r'$T_2^{u,(3)}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$' ,  r'$T_2^{g,(3)}$',
                     r'$E^{u,(1)}$', r'$E^{u,(2)}$' , r'$E^{g,(1)}$', r'$E^{g,(2)}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(3,2,1)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.03)
    sec.set_xticks([0.5, 1.5,3.5,6.5,9.5,12.5,15.5,18.5,20.5,21.5,23.5,26.5,29.5,32.5,35.5,38.5,41,43,45,47], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(28,24, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,2,5,8,11,14,17,20,21,22,25,28,31,34,37,40,42,44,46,47.985], labels = [])
    sec2.tick_params('x',length = 53, width = 1)
    
    # change size of labels on colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 40)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi321_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(3,2,1)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi321_heatmap_yaxis_idx.png")
    
    plt.close()

    ## 3 2 2 ## 
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi322_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi322_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$',                                        
                                        r'$\hat{\epsilon}_1^{\,u}$',r'$\hat{\epsilon}_2^{\,u}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^u$', r'$T_2^{u}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$', r'$E^{u}$',  r'$E^{g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(3,2,2)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.075)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 53, width = 1)

    # change size of labels on colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi322_heatmap.png",dpi = 100)
    
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(3,2,2)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi322_heatmap_yaxis_idx.png")
    plt.close()

    ############################################
    #    different x axis

    #  2 1 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi210_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi210_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,g}$',
                                        r'$\hat{\tau}_1^{\,u,(1)}$',r'$\hat{\tau}_2^{\,u,(1)}$',r'$\hat{\tau}_3^{\,u,(1)}$',
                                        r'$\hat{\tau}_1^{\,u,(2)}$',r'$\hat{\tau}_2^{\,u,(2)}$',r'$\hat{\tau}_3^{\,u,(2)}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
                                        r'$\hat{\epsilon}_1^{\,g,(1)}$',r'$\hat{\epsilon}_2^{\,g,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(2)}$',r'$\hat{\epsilon}_2^{\,g,(2)}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^g$', r'$T_2^{u,(1)}$', r'$T_2^{u,(2)}$', r'$T_2^g$', r'$E^{g,(1)}$',  r'$E^{g,(2)}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(2,1,0)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=0)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23])
    sec.set_xticklabels([l + "    " for l in labels_irreps],rotation = 90,size = 20)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi210_heatmap_xaxis_2.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    # ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_yticklabels(["1","","","","","6","","","","","","12","","","","","","18",
                        "","","","","","24"], size = 15)
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(2,1,0)}$", size = 20)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi210_heatmap_xaxis_2_yaxis_idx.png")
    plt.close()

    ##  2 1 1  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi211_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi211_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$',                                        
                                        r'$\hat{\epsilon}_1^{\,u}$',r'$\hat{\epsilon}_2^{\,u}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^u$', r'$T_2^{u}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$', r'$E^{u}$',  r'$E^{g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(2,1,1)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=0)    
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23])#, labels = labels_irreps)
    sec.set_xticklabels([l + "    " for l in labels_irreps],rotation = 90,size = 20)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi211_heatmap.png",dpi = 100)
    
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(2,1,1)}$",size = 20)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi211_heatmap_xaxis_2_yaxis_idx.png")
    plt.close()

    # 2 2 1 ##

    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi221_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi221_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$',                                        
                                        r'$\hat{\epsilon}_1^{\,u}$',r'$\hat{\epsilon}_2^{\,u}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^u$', r'$T_2^{u}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$', r'$E^{u}$',  r'$E^{g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(2,2,1)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.0)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23])#, labels = labels_irreps)
    sec.set_xticklabels([l + "    " for l in labels_irreps],rotation = 90,size = 20)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi221_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(2,2,1)}$",size = 20)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi221_heatmap_xaxis_2_yaxis_idx.png")
    
    plt.close()
    
    # 3 1 0 ## 

    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi310_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi310_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,g}$',
                                        r'$\hat{\tau}_1^{\,u,(1)}$',r'$\hat{\tau}_2^{\,u,(1)}$',r'$\hat{\tau}_3^{\,u,(1)}$',
                                        r'$\hat{\tau}_1^{\,u,(2)}$',r'$\hat{\tau}_2^{\,u,(2)}$',r'$\hat{\tau}_3^{\,u,(2)}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
                                        r'$\hat{\epsilon}_1^{\,g,(1)}$',r'$\hat{\epsilon}_2^{\,g,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(2)}$',r'$\hat{\epsilon}_2^{\,g,(2)}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^g$', r'$T_2^{u,(1)}$', r'$T_2^{u,(2)}$', r'$T_2^g$', r'$E^{g,(1)}$',  r'$E^{g,(2)}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(3,1,0)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.0)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23])#, labels = labels_irreps)
    sec.set_xticklabels([l + "    " for l in labels_irreps],rotation = 90,size = 20)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi310_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(3,1,0)}$",size = 20)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi310_heatmap_xaxis_2_yaxis_idx.png")
    
    plt.close()

    ##  3 1 1  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi311_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi311_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$',                                        
                                        r'$\hat{\epsilon}_1^{\,u}$',r'$\hat{\epsilon}_2^{\,u}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^u$', r'$T_2^{u}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$', r'$E^{u}$',  r'$E^{g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(3,1,1)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.0)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23])#, labels = labels_irreps)
    sec.set_xticklabels([l + "    " for l in labels_irreps],rotation = 90,size = 20)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi311_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(3,1,1)}$",size = 20)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi311_heatmap_xaxis_2_yaxis_idx.png")
    
    plt.close()

    ##3 2 1 ##

    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi321_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi321_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,u}$',
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,u,(3)}$',r'$\hat{y}^{\,u,(3)}$',r'$\hat{z}^{\,u,(3)}$',
                                        r'$\hat{x}^{\,g,(1)}$',r'$\hat{y}^{\,g,(1)}$',r'$\hat{z}^{\,g,(1)}$',
                                        r'$\hat{x}^{\,g,(2)}$',r'$\hat{y}^{\,g,(2)}$',r'$\hat{z}^{\,g,(2)}$',
                                        r'$\hat{x}^{\,g,(3)}$',r'$\hat{y}^{\,g,(3)}$',r'$\hat{z}^{\,g,(3)}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{a}_2^{\,g}$',
                                        r'$\hat{\tau}_1^{\,u,(1)}$',r'$\hat{\tau}_2^{\,u,(1)}$',r'$\hat{\tau}_3^{\,u,(1)}$',
                                        r'$\hat{\tau}_1^{\,u,(2)}$',r'$\hat{\tau}_2^{\,u,(2)}$',r'$\hat{\tau}_3^{\,u,(2)}$',
                                        r'$\hat{\tau}_1^{\,u,(3)}$',r'$\hat{\tau}_2^{\,u,(3)}$',r'$\hat{\tau}_3^{\,u,(3)}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$', 
                                        r'$\hat{\tau}_1^{\,g,(3)}$',r'$\hat{\tau}_2^{\,g,(3)}$',r'$\hat{\tau}_3^{\,g,(3)}$',                                       
                                        r'$\hat{\epsilon}_1^{\,u,(1)}$',r'$\hat{\epsilon}_2^{\,u,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,u,(2)}$',r'$\hat{\epsilon}_2^{\,u,(2)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(1)}$',r'$\hat{\epsilon}_2^{\,g,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(2)}$',r'$\hat{\epsilon}_2^{\,g,(2)}$'
                                        ]
    labels_irreps = [r'$A_1^u$',r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,  r'$T_1^{u,(3)}$', r'$T_1^{g,(1)}$', r'$T_1^{g,(2)}$' ,  r'$T_1^{g,(3)}$',
                     r'$A_2^u$', r'$A_2^g$', r'$T_2^{u,(1)}$', r'$T_2^{u,(2)}$' ,  r'$T_2^{u,(3)}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$' ,  r'$T_2^{g,(3)}$',
                     r'$E^{u,(1)}$', r'$E^{u,(2)}$' , r'$E^{g,(1)}$', r'$E^{g,(2)}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(3,2,1)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.0)
    special_irrep_labels = [l + "   " for l in labels_irreps]
    special_irrep_labels[0] = labels_irreps[0] + r" $\longrightarrow$  " 
    special_irrep_labels[1] = labels_irreps[1] + r"$\rightarrow$  " 
    special_irrep_labels[8] = labels_irreps[8] + r" $\longrightarrow$  " 
    special_irrep_labels[9] = labels_irreps[9] + r"$\rightarrow$  "
    sec.set_xticks([0.5, 1.5,3.5,6.5,9.5,12.5,15.5,18.5,20.5,21.5,23.5,26.5,29.5,32.5,35.5,38.5,41,43,45,47], 
                   labels = special_irrep_labels, rotation = 90, size = 40)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(28,24, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,2,5,8,11,14,17,20,21,22,25,28,31,34,37,40,42,44,46,47.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)
    
    # change size of labels on colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 40)
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi321_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ylabels = ["" for x in range(len(vecs))]
    ylabels[0] = str(1)
    ylabels[5] = str(6)
    ylabels[11] = str(12)
    ylabels[17] = str(18)
    ylabels[23] = str(24)
    ylabels[29] = str(30)
    ylabels[35] = str(36)
    ylabels[41] = str(42)
    ylabels[47] = str(48)
    ax.set_yticklabels(ylabels,size = 40)
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(3,2,1)}$",size = 40)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi321_heatmap_xaxis_2_yaxis_idx.png")
    
    plt.close()

    ## 3 2 2 ## 
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi322_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopi322_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$',                                        
                                        r'$\hat{\epsilon}_1^{\,u}$',r'$\hat{\epsilon}_2^{\,u}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^u$', r'$T_2^{u}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$', r'$E^{u}$',  r'$E^{g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(3,2,2)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.0)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23])#, labels = labels_irreps)
    sec.set_xticklabels([l + "    " for l in labels_irreps],rotation = 90,size = 20)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)

    # change size of labels on colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi322_heatmap.png",dpi = 100)
    
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(3,2,2)}$",size = 20)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi322_heatmap_xaxis_2_yaxis_idx.png")
    plt.close()

    ###############################################################

    ## p 0 0 ##

    ##  1 0 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopip00_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopip00_vecdata_irreps_seq.npy")
    # print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")

    labels_v = [r'$\hat{a}_1^{\,g}$',r'$\hat{x}^{\,u}$',r'$\hat{y}^{\,u}$',r'$\hat{z}^{\,u}$',r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(p,0,0)}' for x in range(len(vecs))]
    labels_irreps = [r'$A_1^g$',r'$T_1^u$',r'$E^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y],rotation = 0)

    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5, 5], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)

    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,5.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopip00_heatmap.png")
    
    ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(p,0,0)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopip00_heatmap_yaxis_idx.png")

    plt.close()

    #  p p 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopipp0_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopipp0_vecdata_irreps_seq.npy")
    print("irreps pp0:",irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u}$',r'$\hat{y}^{\,u}$',r'$\hat{z}^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(p,p,0)}' for x in range(len(vecs))]
    labels_irreps = [r'$A_1^g$',r'$T_1^u$',r'$T_2^u$', r'$T_2^g$',r'$E^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y], rotation = 0)
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,11], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopipp0_heatmap.png")
    
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(p,p,0)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopipp0_heatmap_yaxis_idx.png")
    plt.close()

    #  p p p  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopippp_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopippp_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u}$',r'$\hat{y}^{\,u}$',r'$\hat{z}^{\,u}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
    ]
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(p,p,p)}' for x in range(len(vecs))]
    labels_irreps = [r'$A_1^g$',r'$T_1^u$',r'$A_2^u$', r'$T_2^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y], rotation = 0)
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5, 4.5, 6.5], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,5,7.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopippp_heatmap.png")
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(p,p,p)}$")
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopippp_heatmap_yaxis_idx.png")
    
    plt.close()

    #  p q 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopipq0_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopipq0_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,g}$',
                                        r'$\hat{\tau}_1^{\,u,(1)}$',r'$\hat{\tau}_2^{\,u,(1)}$',r'$\hat{\tau}_3^{\,u,(1)}$',
                                        r'$\hat{\tau}_1^{\,u,(2)}$',r'$\hat{\tau}_2^{\,u,(2)}$',r'$\hat{\tau}_3^{\,u,(2)}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
                                        r'$\hat{\epsilon}_1^{\,g,(1)}$',r'$\hat{\epsilon}_2^{\,g,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(2)}$',r'$\hat{\epsilon}_2^{\,g,(2)}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^g$', r'$T_2^{u,(1)}$', r'$T_2^{u,(2)}$', r'$T_2^g$', r'$E^{g,(1)}$',  r'$E^{g,(2)}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(p,q,0)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=0)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23])
    sec.set_xticklabels([l + "    " for l in labels_irreps],rotation = 90,size = 20)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi210_heatmap_xaxis_2.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    # ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_yticklabels(["1","","","","","6","","","","","","12","","","","","","18",
                        "","","","","","24"], size = 15)
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(q,p,0)}$", size = 20)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopipq0_heatmap_xaxis_2_yaxis_idx.png")
    plt.close()

    ##  q p p  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopippq_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopippq_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{\tau}_1^{\,u}$',r'$\hat{\tau}_2^{\,u}$',r'$\hat{\tau}_3^{\,u}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$',                                        
                                        r'$\hat{\epsilon}_1^{\,u}$',r'$\hat{\epsilon}_2^{\,u}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,r'$T_1^{g}$',r'$A_2^u$', r'$T_2^{u}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$', r'$E^{u}$',  r'$E^{g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(q,p,p)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=0)    
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,10.5, 12.5, 15.5,18.5,21,23])#, labels = labels_irreps)
    sec.set_xticklabels([l + "    " for l in labels_irreps],rotation = 90,size = 20)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi211_heatmap.png",dpi = 100)
    
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(q,p,p)}$",size = 20)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopippq_heatmap_xaxis_2_yaxis_idx.png")
    plt.close()

    ##r q p ##

    vecs = np.load("D:/Master/Masterarbeit/results/twopi/data/twopipqr_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi/data/twopipqr_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,u}$',
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,u,(1)}$',r'$\hat{y}^{\,u,(1)}$',r'$\hat{z}^{\,u,(1)}$',
                                        r'$\hat{x}^{\,u,(2)}$',r'$\hat{y}^{\,u,(2)}$',r'$\hat{z}^{\,u,(2)}$',
                                        r'$\hat{x}^{\,u,(3)}$',r'$\hat{y}^{\,u,(3)}$',r'$\hat{z}^{\,u,(3)}$',
                                        r'$\hat{x}^{\,g,(1)}$',r'$\hat{y}^{\,g,(1)}$',r'$\hat{z}^{\,g,(1)}$',
                                        r'$\hat{x}^{\,g,(2)}$',r'$\hat{y}^{\,g,(2)}$',r'$\hat{z}^{\,g,(2)}$',
                                        r'$\hat{x}^{\,g,(3)}$',r'$\hat{y}^{\,g,(3)}$',r'$\hat{z}^{\,g,(3)}$',
                                        r'$\hat{a}_2^{\,u}$',
                                        r'$\hat{a}_2^{\,g}$',
                                        r'$\hat{\tau}_1^{\,u,(1)}$',r'$\hat{\tau}_2^{\,u,(1)}$',r'$\hat{\tau}_3^{\,u,(1)}$',
                                        r'$\hat{\tau}_1^{\,u,(2)}$',r'$\hat{\tau}_2^{\,u,(2)}$',r'$\hat{\tau}_3^{\,u,(2)}$',
                                        r'$\hat{\tau}_1^{\,u,(3)}$',r'$\hat{\tau}_2^{\,u,(3)}$',r'$\hat{\tau}_3^{\,u,(3)}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$', 
                                        r'$\hat{\tau}_1^{\,g,(3)}$',r'$\hat{\tau}_2^{\,g,(3)}$',r'$\hat{\tau}_3^{\,g,(3)}$',                                       
                                        r'$\hat{\epsilon}_1^{\,u,(1)}$',r'$\hat{\epsilon}_2^{\,u,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,u,(2)}$',r'$\hat{\epsilon}_2^{\,u,(2)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(1)}$',r'$\hat{\epsilon}_2^{\,g,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(2)}$',r'$\hat{\epsilon}_2^{\,g,(2)}$'
                                        ]
    labels_irreps = [r'$A_1^u$',r'$A_1^g$',r'$T_1^{u,(1)}$', r'$T_1^{u,(2)}$' ,  r'$T_1^{u,(3)}$', r'$T_1^{g,(1)}$', r'$T_1^{g,(2)}$' ,  r'$T_1^{g,(3)}$',
                     r'$A_2^u$', r'$A_2^g$', r'$T_2^{u,(1)}$', r'$T_2^{u,(2)}$' ,  r'$T_2^{u,(3)}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$' ,  r'$T_2^{g,(3)}$',
                     r'$E^{u,(1)}$', r'$E^{u,(2)}$' , r'$E^{g,(1)}$', r'$E^{g,(2)}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(r,q,p)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.0)
    special_irrep_labels = [l + "   " for l in labels_irreps]
    special_irrep_labels[0] = labels_irreps[0] + r" $\longrightarrow$  " 
    special_irrep_labels[1] = labels_irreps[1] + r"$\rightarrow$  " 
    special_irrep_labels[8] = labels_irreps[8] + r" $\longrightarrow$  " 
    special_irrep_labels[9] = labels_irreps[9] + r"$\rightarrow$  "
    sec.set_xticks([0.5, 1.5,3.5,6.5,9.5,12.5,15.5,18.5,20.5,21.5,23.5,26.5,29.5,32.5,35.5,38.5,41,43,45,47], 
                   labels = special_irrep_labels, rotation = 90, size = 40)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(28,24, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,2,5,8,11,14,17,20,21,22,25,28,31,34,37,40,42,44,46,47.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)
    
    # change size of labels on colorbar
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 40)
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi321_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ylabels = ["" for x in range(len(vecs))]
    ylabels[0] = str(1)
    ylabels[5] = str(6)
    ylabels[11] = str(12)
    ylabels[17] = str(18)
    ylabels[23] = str(24)
    ylabels[29] = str(30)
    ylabels[35] = str(36)
    ylabels[41] = str(42)
    ylabels[47] = str(48)
    ax.set_yticklabels(ylabels,size = 40)
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(r,q,q)}$",size = 40)
    fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopipqr_heatmap_xaxis_2_yaxis_idx.png")
    
    plt.close()