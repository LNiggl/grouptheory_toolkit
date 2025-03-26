import numpy as np
import seaborn as sns
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
    

    ## p 0 0 ##

    vecs = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_p00_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_p00_vecdata_irreps_seq.npy")
    # print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()

    labels_v = [r'$\hat{a}_1^{\,g}$',r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(p,0,0)}' for x in range(len(vecs))]
    labels_irreps = [r'$A_1^g$',r'$E^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y],rotation = 0)

    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)

    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,2.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)
    # plt.show()
    
    # fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist_p00_heatmap.png",dpi = 400)
    
    ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(p,0,0)}$")
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist_p00_heatmap_yaxis_idx.png",dpi = 400)

    plt.close()

    #  p p 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_pp0_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_pp0_vecdata_irreps_seq.npy")
    print("irreps pp0:",irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(p,p,0)}' for x in range(len(vecs))]
    labels_irreps = [r'$A_1^g$',r'$T_2^g$',r'$E^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y], rotation = 0)
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5, 5], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,5.985],labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    # plt.show()
    # fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist_pp0_heatmap.png",dpi = 400)
    
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(p,p,0)}$")
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist_pp0_heatmap_yaxis_idx.png",dpi = 400)
    plt.close()

    ## p p p  ##

    vecs = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_ppp_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_ppp_vecdata_irreps_seq.npy")
    # print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()

    labels_v = [r'$\hat{a}_1^{\,g}$',r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(p,p,p)}' for x in range(len(vecs))]
    labels_irreps = [r'$A_1^g$',r'$T_2^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y],rotation = 0)

    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)

    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,5.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    
    # fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist_ppp_heatmap.png",dpi = 400)
    
    ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(p,p,p)}$")
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist_ppp_heatmap_yaxis_idx.png",dpi = 400)

    plt.close()


    #  p q 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_pq0_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_pq0_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{a}_2^{\,g}$',
                                        r'$\hat{\tau}_1^{\,g}$',r'$\hat{\tau}_2^{\,g}$',r'$\hat{\tau}_3^{\,g}$',
                                        r'$\hat{\epsilon}_1^{\,g,(1)}$',r'$\hat{\epsilon}_2^{\,g,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(2)}$',r'$\hat{\epsilon}_2^{\,g,(2)}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{g}$',r'$A_2^g$',  r'$T_2^g$', r'$E^{g,(1)}$',  r'$E^{g,(2)}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(p,q,0)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()


    sec = ax.secondary_xaxis(location=0)
    sec.set_xticks([0.5, 2.5, 4.5, 6.5,9,11])
    sec.set_xticklabels(["\n" + l for l in labels_irreps],size = 20)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,5,8,10,11.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)

    # sec = ax.secondary_xaxis(location=-0.1)
    # sec.set_xticks([0.5, 2.5, 4.5, 6.5,9,11], labels = labels_irreps)
    # sec.tick_params('x',length = 0)
    # fig.subplots_adjust(bottom = 0.2,left = 0.2)

    # sec2 = ax.secondary_xaxis(location = 0)
    # sec2.set_xticks([0.015,1,4,5,8,10,11.985], labels = [])
    # sec2.tick_params('x',length = 40, width = 1)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    
    # fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist210_heatmap_xaxis_2.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))],rotation = 0, size = 15)
    # ax.set_yticklabels(["1","","","","","6","","","","","","12","","","","","","18",
                        # "","","","","","24"], size = 15)
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(q,p,0)}$", size = 20)
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist_pq0_heatmap_xaxis_2_yaxis_idx.png",dpi = 400)
    plt.close()

    ##  q p p  ##
    vecs = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_ppq_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_ppq_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        r'$\hat{a}_1^{\,g}$',
                                        r'$\hat{x}^{\,g}$',r'$\hat{y}^{\,g}$',r'$\hat{z}^{\,g}$',
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$',                                        
                                        r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{g}$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$',  r'$E^{g}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(q,p,p)}' for x in range(len(vecs))]
    # ax.tick_params(axis = 'both', which = 'major', labelsize = 10)
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)
    ax.set_yticks([x + 0.5 for x in range(len(vecs))], labels = [r'${}$'.format(l) for l in labels_y])
    # plt.show()
    sec = ax.secondary_xaxis(location=0)    
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,11], labels = [])#, 12.5, 15.5,18.5,21,23])#, labels = labels_irreps)
    sec.set_xticklabels(["\n" + l for l in labels_irreps],rotation = 0,size = 20)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    fig.set_size_inches(12,10, forward= True)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11.985],labels = [])#,14,17,20,22,23.985], labels = [])
    sec2.tick_params('x',length = 15, width = 2)

    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize = 15)
    # plt.show()
    # fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist211_heatmap.png",dpi = 100)
    
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    ax.set_yticklabels([str(x+1) for x in range(len(vecs))], rotation = 0, size = 15)
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(q,p,p)}$",size = 20)
    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist_ppq_heatmap_xaxis_2_yaxis_idx.png",dpi = 400)
    plt.close()

    ##r q p ##

    vecs = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_pqr_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/twopi_indist/data/twopi_indist_pqr_vecdata_irreps_seq.npy")
    print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")
    # plt.show()
    labels_v = [
                                        
                                        r'$\hat{a}_1^{\,g}$',
                                        
                                        r'$\hat{x}^{\,g,(1)}$',r'$\hat{y}^{\,g,(1)}$',r'$\hat{z}^{\,g,(1)}$',
                                        r'$\hat{x}^{\,g,(2)}$',r'$\hat{y}^{\,g,(2)}$',r'$\hat{z}^{\,g,(2)}$',
                                        r'$\hat{x}^{\,g,(3)}$',r'$\hat{y}^{\,g,(3)}$',r'$\hat{z}^{\,g,(3)}$',
                                        
                                        r'$\hat{a}_2^{\,g}$',
                                        
                                        r'$\hat{\tau}_1^{\,g,(1)}$',r'$\hat{\tau}_2^{\,g,(1)}$',r'$\hat{\tau}_3^{\,g,(1)}$',
                                        r'$\hat{\tau}_1^{\,g,(2)}$',r'$\hat{\tau}_2^{\,g,(2)}$',r'$\hat{\tau}_3^{\,g,(2)}$', 
                                        r'$\hat{\tau}_1^{\,g,(3)}$',r'$\hat{\tau}_2^{\,g,(3)}$',r'$\hat{\tau}_3^{\,g,(3)}$',                                       
                                       
                                        r'$\hat{\epsilon}_1^{\,g,(1)}$',r'$\hat{\epsilon}_2^{\,g,(1)}$',
                                        r'$\hat{\epsilon}_1^{\,g,(2)}$',r'$\hat{\epsilon}_2^{\,g,(2)}$'
                                        ]
    labels_irreps = [r'$A_1^g$',r'$T_1^{g,(1)}$', r'$T_1^{g,(2)}$' ,  r'$T_1^{g,(3)}$',
                    r'$A_2^g$', r'$T_2^{g,(1)}$', r'$T_2^{g,(2)}$' ,  r'$T_2^{g,(3)}$',
                    r'$E^{g,(1)}$', r'$E^{g,(2)}$']
    labels_y = [r'\hat{b}' + '_{'+ f'{x+1}' +'}^{(r,q,p)}' for x in range(len(vecs))]
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = [])#labels_v)

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
    
    # fig.savefig("D:/Master/Masterarbeit/results/twopi/plots/twopi210_heatmap_xaxis_2.png",dpi = 100)
     ## yaxis by indices ##
    ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    # ax.set_yticklabels([str(x+1) for x in range(len(vecs))])
    ax.set_yticklabels(["1","","","","","6","","","","","","12","","","","","","18",
                        "","","","","","24"], size = 15)
    ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(r,q,p)}$", size = 20)
    plt.show()
    
    # fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist_pqr_heatmap.png",dpi = 100)
     ## yaxis by indices ##
    # ax.set_yticks([x + 0.5 for x in range(len(vecs))])
    # ylabels = ["" for x in range(len(vecs))]
    # ylabels[0] = str(1)
    # ylabels[5] = str(6)
    # ylabels[11] = str(12)
    # ylabels[17] = str(18)
    # ylabels[23] = str(24)
    # ylabels[29] = str(30)
    # ylabels[35] = str(36)
    # ylabels[41] = str(42)
    # ylabels[47] = str(48)
    # ax.set_yticklabels(ylabels,size = 40)
    # ax.set_ylabel(r"index $i$ : $\hat{b}_i^{(r,q,p)}$",size = 40)
    fig.savefig("D:/Master/Masterarbeit/results/twopi_indist/plots/twopi_indist_pqr_heatmap_xaxis_2_yaxis_idx.png",dpi = 400)
    
    plt.close()