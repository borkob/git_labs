The minimization of the autocorrelation (energy) function, or equivalently the maximization of the binary 
sequence merit factor, has important applications in communication engineering. To physicists, the optimum
solution of the labs problem corresponds to the ground state of a generalized one-dimensional Ising spin system
with long range 4-spin interactions, also known as the Bernasconi model with **aperiodic autocorrelation**. 
Finding optimum sequences of the labs problem is significantly harder than solving the Ising spin-glass problems
with limited interaction and periodic boundary conditions. 

In fact, the runtime asymptotic complexity of the labs problem is much harder than the Optimal Golomb Ruler problem
that continues to be solved under massively parallel computation for the past few decades 
(see http://en.wikipedia.org/wiki/Golomb_ruler and http://www.distributed.net/OGR).

A crowd-sourcing server prototype to facilitate experimentation and push the frontiers on finding new **best-known values** (BKVs), for the labs problem is under construction. Researchers will be able to download the labs solver lssOrelE with the data set and simply **run the solver on their host** in the browser. Upon finding the solution (or a breakthrough value that improves the current **BKV**), the browser will return solution details to the server and start a new run on the local host either with the current or the new **BKV**. For an interactive labs puzzle and the timeline of the crowd-sourcing server prototype availability see http://labraj.feri.um.si/en/B.labs.

Please use this citation when referring to this archive:

B. Bošković, F. Brglez and J. Brest, **A GitHub Archive for Solvers and Solutions of the labs problem**, 2016


@misc{OPUS2-git_labs-Boskovic,

	Author = {B. Bo\v{s}kovi\'{c} and F. Brglez and J. Brest},

	Howpublished = {{For updates, see \url{https://github.com/borkob/git_labs}. }},

	Month = {January},

	Title = {{A GitHub Archive for Solvers and Solutions of the {\tt labs} problem}},

	Year = {2016}

}

The organization of this archive, including the open-source state-of-the-art solvers and summaries of experimental results,
is based on this paper:


B. Bošković, F. Brglez and J. Brest, **Low-Autocorrelation Binary Sequences: On Improved Merit Factors and Runtime Predictions to Achieve Them**, http://arxiv.org/abs/1406.5301, under review

@article{OPUS2-labs-2016-arxiv-Boskovic,

	Author = {B. Bo\v{s}kovi\'{c} and F. Brglez and J. Brest},

	Journal = {http://arxiv.org/abs/1406.5301, also under journal review}},

	Title = {{Low-Autocorrelation Binary Sequences: On Improved Merit Factors and Runtime Predictions to Achieve Them}},

	Year = {2016}

}

