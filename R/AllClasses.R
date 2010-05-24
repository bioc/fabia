#
#
# Author: SEPP HOCHREITER
###############################################################################


setClass('Factorization',
         representation = representation(
           parameters = 'list',
           n = 'numeric',
           p1 = 'numeric',
           p2 = 'numeric',
           l = 'numeric',
           center = 'vector',
           scaleData = 'vector',
           X = 'matrix',
           L = 'matrix',
           Z = 'matrix',
           M = 'matrix',
           LZ = 'matrix',
           U = 'matrix',
           avini = 'vector',
           xavini = 'vector',
           ini = 'matrix',
           Psi = 'vector',
           lapla = 'matrix'))
