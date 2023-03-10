The git repository (available at https://github.com/charliesire/quantization_Campbell2D.git) provides the codes and data to reproduce all the experiments related to the Campbell2D function that are described in the article. 
More precisely : 
 - FunQuant-0.1.3 is the version 0.1.3 of the package FunQuant, that contains the functions to perform the Prototype Maps Algorithm, to tune the metamodel and compute the performance metrics.
 - GpOutput2D-main contains the code from Elodie Perrin to perform FPCA combined with Gaussian Processes. 
 - Campbell2D.R is the Campbell2D function generating the Campbell maps.
 - NewFitting_Charlie_v090821.RData are historical data related to the offshore conditions, providing the probabilistic distributions.
 - Campbell_utils.R contains different functions useful for all the notebooks.
 - PMalgo_true.Rmd performs the lloyd algorithm with the true campbell maps.
 - perf_probas.Rmd tunes the hyperparameters of the metamodel and evaluates the relative probability error.
 - PMalgo_pred.Rmd performs the lloyd algorithm with the predicted campbell maps.
 - compute_probas.Rmd computes the probabilities associated to the Voronoi cells with a higher number of maps that the one used in the lloyd algorithm to increase precision.
 - importance_sampling_error.Rmd compute the errors from the importance sampling :
    +	The IS coefficient of variation of the membership probability
    + The IS centroid standard deviation
  - error_quanti.Rmd

It also contains a file is_perf_metric.pdf providing details about the computation of the importance sampling performance metrics, and a file complexity_pma.pdf explaining the calculation of the complexity of the Prototype Maps Algorithm.
