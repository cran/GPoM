#' @title GPoM package: Generalized Polynomial Modelling
#' @description GPoM is a platform dedicated to the Global Modelling
#' technique wjich aim is to obtain deterministic models of
#' Ordinary Differential Equations from observational time series.
#' It can be used to detect couplings between observed variables
#' and for the inference of causal networks.
#' It can also be applied to single time series.
#' The present package focuses on models in Ordinary Differential
#' Equations of polynomial form.
#' Although designed for modelling chaos, it can be also applied
#' to any linear or nonlinear dynamical system for which
#' a deterministic behavior is expected.
#'
#' @author Sylvain Mangiarotti, Flavie Le Jean,
#' Malika Chassan, Laurent Drapeau, Mireille Huc.
#'
#' Maintainer: M. Huc <mireille.huc@@cesbio.cnes.fr>
#'
#' @references
#' [1] J. P. Crutchfield and B. S. McNamara, 1987.
#' Equations of motion from a data series,
#' Complex Syst. 1, 417-452. \cr
#' [2] Gouesbet G., Letellier C., 1994.
#' Global vector-field reconstruction by using a multivariate
#' polynomial L2 approximation on nets,
#' Physical Review E, 49 (6), 4955-4972. \cr
#' [3] C. Letellier, L. Le Sceller, E. Marechal, P. Dutertre, B. Maheu,
#' G. Gouesbet, Z. Fei, and J. L. Hudson, 1995.
#' Global vector field reconstruction from a chaotic experimental
#' signal in copper electrodissolution,
#' Phys. Rev. E 51, 4262-4266. \cr
#' [4] J. Maquet, C. Letellier, and L. A. Aguirre, 2007.
#' Global models from the Canadian Lynx cycles as a first evidence
#' for chaos in real ecosystems,
#' J. Math. Biol. 55(1), 21-39. \cr
#' [5] Mangiarotti S., Coudret R., Drapeau L., & Jarlan L., 2012.
#' Polynomial search and global modeling : Two algorithms for
#' modeling chaos,
#' Physical Review E, 86, 046205. \cr
#' [6] Mangiarotti S., Drapeau L. & Letellier C., 2014.
#' Two chaotic models for cereal crops observed from satellite in
#' northern Morocco.
#' Chaos, 24(2), 023130. \cr
#' [7] Mangiarotti S., 2015. Low dimensional chaotic models for the
#' plague epidemic in Bombay (1896-1911).
#' Chaos, Solitons and Fractals, 81A, 184-186. \cr
#' [8] Mangiarotti S., Peyre M. & Huc M.,
#' A chaotic model for the epidemic of Ebola Virus Disease in West Africa (2013-2016).
#' Chaos, 26, 113112, 2016. \cr
#' [9] Mangiarotti S., 2014. Modelisation globale et Caracterisation
#' Topologique de dynamiques environnementales - de l'analyse des
#' enveloppes fluides et du couvert de surface de la Terre a la
#' caracterisation topolodynamique du chaos.
#' Habilitation to Direct Research,
#' University of Toulouse 3, France. \cr
#' @note
#' FOR USERS \cr
#' This package was developped in Centre d'Etudes Spatiales de
#' la Biosphere http://www.cesbio.ups-tlse.fr/.
#' If you want to apply this package to single time series,
#' please quote [5].
#' If you want to apply it to multivariate time series,
#' please quote [9].
#' If you want to apply it to infer causal relations among
#' time series, please quote [7].
#' @note
#' HISTORICAL BACKGROUND \cr
#' The global modelling technique was initiated during
#' the early 1990s [1-3]. It takes its background from
#' the Theory of Nonlinear Dynamical Systems.
#' The approach became applicable to the analysis of real world
#' environmental behaviours by the end of the 2000s [4-6].
#' Recent works have shown that the approach could be applied to
#' numerous other dynamical behaviors [7-9].
#' Global modelling aims to obtain deterministic models
#' directly from observed time series.
#'
#' @docType package
#' @name GPoM
#' @keywords global modeling, nonlinear dynamical systems, chaos,
#' causal inference, time series analysis, data learning
NA