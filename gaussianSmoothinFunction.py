def gaussianSmoothing2D(x_vals,y_vals,sigma):
    y_smth = np.zeros(y_vals.shape) 
    for it in range(0,len(x_vals)):
        x_position      = x_vals[it]
        gaussian_kernel = np.exp(-(x_vals - x_position) ** 2 / (2 * sigma ** 2))
        gaussian_kernel = gaussian_kernel / np.sum(gaussian_kernel)
        y_smth[it]      = np.sum(y_vals * gaussian_kernel)
    return(y_smth)
