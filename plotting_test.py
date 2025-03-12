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
    vecs = np.load("D:/Master/Masterarbeit/results/test_01/twopi100_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/test_01/twopi100_vecdata_irreps_seq.npy")
    # print(irreps)

    fig,ax = plt.subplots()
    ax = sns.heatmap(vecs.T,linewidth = 0.5, cmap = "coolwarm")

    labels_v = [r'$\hat{a}_1^{\,g}$',r'$\hat{x}^{\,u}$',r'$\hat{y}^{\,u}$',r'$\hat{z}^{\,u}$',r'$\hat{\epsilon}_1^{\,g}$',r'$\hat{\epsilon}_2^{\,g}$']
    labels_irreps = [r'$A_1^g$',r'$T_1^u$',r'$E^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)

    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5, 5], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)

    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,5.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/test_01/twopi100_heatmap.png")
    plt.close()

    #  1 1 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/test_01/twopi110_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/test_01/twopi110_vecdata_irreps_seq.npy")
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
    labels_irreps = [r'$A_1^g$',r'$T_1^u$',r'$T_2^u$', r'$T_2^g$',r'$E^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5, 5.5, 8.5,11], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,7,10,11.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/test_01/twopi110_heatmap.png")
    plt.close()

    #  1 1 1  ##
    vecs = np.load("D:/Master/Masterarbeit/results/test_01/twopi111_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/test_01/twopi111_vecdata_irreps_seq.npy")
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
    labels_irreps = [r'$A_1^g$',r'$T_1^u$',r'$A_2^u$', r'$T_2^g$']
    ax.set_xticks([x + 0.5 for x in range(len(vecs))], labels = labels_v)
    # plt.show()
    sec = ax.secondary_xaxis(location=-0.1)
    sec.set_xticks([0.5, 2.5, 4.5, 6.5], labels = labels_irreps)
    sec.tick_params('x',length = 0)
    fig.subplots_adjust(bottom = 0.2,left = 0.2)
    # plt.show()
    sec2 = ax.secondary_xaxis(location = 0)
    sec2.set_xticks([0.015,1,4,5,7.985], labels = [])
    sec2.tick_params('x',length = 40, width = 1)

    plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/test_01/twopi111_heatmap.png")
    plt.close()

    #  2 1 0  ##
    vecs = np.load("D:/Master/Masterarbeit/results/test_01/twopi210_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/test_01/twopi210_vecdata_irreps_seq.npy")
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

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/test_01/twopi210_heatmap.png",dpi = 100)
    plt.close()

    ##  2 1 1  ##
    vecs = np.load("D:/Master/Masterarbeit/results/test_01/twopi211_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/test_01/twopi211_vecdata_irreps_seq.npy")
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

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/test_01/twopi211_heatmap.png",dpi = 100)
    plt.close()

    # 2 2 1 ##

    vecs = np.load("D:/Master/Masterarbeit/results/test_01/twopi221_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/test_01/twopi221_vecdata_irreps_seq.npy")
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

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/test_01/twopi221_heatmap.png",dpi = 100)
    plt.close()
    
    # 3 1 0 ## 

    vecs = np.load("D:/Master/Masterarbeit/results/test_01/twopi310_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/test_01/twopi310_vecdata_irreps_seq.npy")
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

    plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/test_01/twopi310_heatmap.png",dpi = 100)
    plt.close()

    ##  3 1 1  ##
    vecs = np.load("D:/Master/Masterarbeit/results/test_01/twopi311_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/test_01/twopi311_vecdata_irreps_seq.npy")
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

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/test_01/twopi311_heatmap.png",dpi = 100)
    plt.close()

    # 3 2 1 ## 
    vecs = np.load("D:/Master/Masterarbeit/results/test_01/twopi321_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/test_01/twopi321_vecdata_irreps_seq.npy")
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

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/test_01/twopi321_heatmap.png",dpi = 100)
    plt.close()

    ## 3 2 2 ## 
    vecs = np.load("D:/Master/Masterarbeit/results/test_01/twopi322_vecdata.npy")
    irreps = np.load("D:/Master/Masterarbeit/results/test_01/twopi322_vecdata_irreps_seq.npy")
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

    # plt.show()
    fig.savefig("D:/Master/Masterarbeit/results/test_01/twopi322_heatmap.png",dpi = 100)
    plt.close()

    ## thoughts: 
    # could adjust such that no line between irrep and basis vectors
    # then omit basis vectors and axes look sort of consistent
    # or adjust size of labels such that letters fit

    #   WHATEVER I CHOOSE: DONT OVERWRITE EXISTING PLOTS!
