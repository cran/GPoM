#' @title GPoM package: Generalized Polynomial Modelling
#' @description GPoM is a platform dedicated to the Global Modelling
#' technique. Its aim is to obtain deterministic models of
#' Ordinary Differential Equations from observational time series.
#' It applies to single and to multiple time series.
#' With single time series, it can be used:
#' to detect low-dimnesional determinism and low-dimensional
#' (deterministic) chaos. It can also be used to
#' characterize the observed behavior, using the obtained
#' models as a proxy of the original dynamics,
#' as far as the model validation could be checked.
#' With multiple time series, it can be used:
#' to detect couplings between observed variables,
#' to infer causal networks,
#' and to reformulate the original equations of the observed
#' system (retro-modelling).
#' The present package focuses on models in Ordinary Differential
#' Equations of polynomial form.
#' The package was designed to model weakly predictable dynamical
#' behaviors (such as chaotic behaviors).
#' Of course, it can also apply to more of fully predictable
#' behavior, either linear or nonlinear.
#' Several vignettes are associated to the package
#' which can be used as a tutorial,
#' and it also provides an overlook of the diversity
#' of applications and at the performances of the tools.
#' Users are kindly asked to quote the corresponding
#' references when using the package (see hereafter).
#'
#' @author Sylvain Mangiarotti, Flavie Le Jean,
#' Malika Chassan, Laurent Drapeau, Mireille Huc.
#'
#' Maintainer: M. Huc <mireille.huc@@u-paris2.fr>
#'
#' @references
#' [1] J. P. Crutchfield and B. S. McNamara, 1987.
#' Equations of motion from a data series,
#' Complex Systems. 1, 417-452. \cr
#' [2] Gouesbet G., Letellier C., 1994.
#' Global vector-field reconstruction by using a multivariate
#' polynomial L2 approximation on nets,
#' Physical Review E, 49 (6), 4955-4972. \cr
#' [3] C. Letellier, L. Le Sceller, E. Marechal, P. Dutertre, B. Maheu,
#' G. Gouesbet, Z. Fei, and J. L. Hudson, 1995.
#' Global vector field reconstruction from a chaotic experimental
#' signal in copper electrodissolution,
#' Physical Review E, 51, 4262-4266. \cr
#' [4] L. A. Aguirre & C. Letellier,
#' Modeling nonlinear dynamics and chaos: A review,
#' Mathematical Problems in Engineering, 2009, 238960. \cr
#' C. Letellier, L. Le Sceller, E. Marechal, P. Dutertre, B. Maheu,
#' G. Gouesbet, Z. Fei, and J. L. Hudson, 1995.
#' Global vector field reconstruction from a chaotic experimental
#' signal in copper electrodissolution,
#' Physical Review E 51, 4262-4266. \cr
#' [5] J. Maquet, C. Letellier, and L. A. Aguirre, 2007.
#' Global models from the Canadian Lynx cycles as a first evidence
#' for chaos in real ecosystems,
#' Juornal of Mathematical Biology. 55(1), 21-39. \cr
#' [6] Mangiarotti S., Coudret R., Drapeau L., & Jarlan L., 2012.
#' Polynomial search and global modeling : Two algorithms for
#' modeling chaos,
#' Physical Review E, 86, 046205. \cr
#' [7] Mangiarotti S., Drapeau L. & Letellier C., 2014.
#' Two chaotic models for cereal crops observed from satellite in
#' northern Morocco.
#' Chaos, 24(2), 023130. \cr
#' [8] Mangiarotti S., 2015. Low dimensional chaotic models for the
#' plague epidemic in Bombay (1896-1911).
#' Chaos, Solitons and Fractals, 81A, 184-186. \cr
#' [9] Mangiarotti S., Peyre M. & Huc M.,
#' A chaotic model for the epidemic of Ebola Virus Disease in West Africa (2013-2016).
#' Chaos, 26, 113112, 2016. \cr
#' [10] Mangiarotti S., 2014. Modelisation globale et Caracterisation
#' Topologique de dynamiques environnementales - de l'analyse des
#' enveloppes fluides et du couvert de surface de la Terre a la
#' caracterisation topolodynamique du chaos.
#' Habilitation to Direct Research,
#' University of Toulouse 3, France. \cr
#' [11] Mangiarotti S., Sharma A.K., Corgne S., Hubert-Moy L., Ruiz L., Sekhar M., Kerr Y.,
#' Can the global modelling technique be used for crop classification?
#' Chaos, Solitons & Fractals, in press. \cr
#'
#' @note
#' FOR USERS \cr
#' This package was developped at Centre d'Etudes Spatiales de
#' la Biosphere (Cesbio, UMR 5126, UPS-CNRS-CNES-IRD,
#' http://www.cesbio.ups-tlse.fr).
#' An important part of the developments were funded by
#' the French program Les Enveloppes Fluides et l'Environnement
#' (LEFE, MANU, projets GloMo, SpatioGloMo and MoMu).
#' The French program Défi InFiNiTi (CNRS) and PNTS
#' are also acknowledged (projects Crops'IChaos and Musc & SlowFast).
#'
#' If you apply this package to single time series,
#' please quote [6].
#' If you apply it to multivariate time series,
#' please quote [10].
#' If you apply it to infer couplings among
#' time series, please quote [8].
#' If you apply it to classification, please quote [11].
#'
#' @note
#' HISTORICAL BACKGROUND \cr
#' The global modelling technique was initiated during
#' the early 1990s [1-3]. It takes its background from
#' the Theory of Nonlinear Dynamical Systems.
#' Earlier investigations can also be found in the fields
#' of Electrical Engineering and Statistics but these
#' mainly focused on linear problems [4].
#' The approach became applicable to the analysis of real world
#' environmental behaviours by the end of the 2000s [5-7].
#' Recent works have shown that the approach could be applied to
#' numerous other dynamical behaviors [8-10].
#' Global modelling aims to obtain deterministic models
#' directly from observed time series.
#'
#' @encoding UTF-8
#'
#' @docType package
#' @name GPoM-package
#' @concept global modeling
#' @concept nonlinear dynamical systems
#' @concept chaos
#' @concept causal inference
#' @concept time series analysis
#' @concept data learning
NA
