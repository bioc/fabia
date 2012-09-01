#
#
# Author: SEPP HOCHREITER
###############################################################################


setClass("Factorization",
         representation = representation(
           parameters = "list",
           n = "numeric",
           p1 = "numeric",
           p2 = "numeric",
           l = "numeric",
           center = "vector",
           scaleData = "vector",
           X = "matrix",
           L = "matrix",
           Z = "matrix",
           M = "matrix",
           LZ = "matrix",
           U = "matrix",
           avini = "vector",
           xavini = "vector",
           ini = "matrix",
           Psi = "vector",
           lapla = "matrix"))


setValidity("Factorization",
    function(object)
    {
        if (!is.list(slot(object, "parameters")))
        {
            return("slot >parameters< must be a list!")
        }
        else if (!is.numeric(slot(object, "n")) || length(slot(object, "n")) != 1 ||
                 slot(object, "n") < 2)
        {
            return("slot >n< must be an integer number larger 1!")
        }
        else if (!is.numeric(slot(object, "p1")) || length(slot(object, "p1")) != 1 ||
                 slot(object, "p1") < 1)
        {
            return("slot >p1< must be an integer number larger 0!")
        }
        else if (!is.numeric(slot(object, "p2")) || length(slot(object, "p2")) != 1 ||
                 slot(object, "p2") < 1)
        {
            return("slot >p2< must be an integer number larger 01")
        }
        else if (!is.numeric(slot(object, "l")) || length(slot(object, "l")) != 1 ||
                 slot(object, "l") < 2)
        {
            return("slot >l< must be an integer number larger 1!")
        }
        else if (!is.vector(slot(object, "center")))
        {
            return("slot >center< must be a vector!")
        }
        else if (!is.vector(slot(object, "scaleData")))
        {
            return("slot >scaleData< must be a vector!")
        }
        else if (!is.matrix(slot(object, "X")))
        {
            return("slot >X< must be a matrix!")
        }
        else if (!is.matrix(slot(object, "L")))
        {
            return("slot >L< must be a matrix!")
        }
        else if (!is.matrix(slot(object, "Z")))
        {
            return("slot >Z< must be a matrix!")
        }
        else if (!is.matrix(slot(object, "M")))
        {
            return("slot >M< must be a matrix!")
        }
        else if (!is.matrix(slot(object, "LZ")))
        {
            return("slot >LZ< must be a matrix!")
        }
        else if (!is.matrix(slot(object, "U")))
        {
            return("slot >U< must be a matrix!")
        }
        else if (!is.vector(slot(object, "avini")))
        {
            return("slot >avini< must be a vector!")
        }
        else if (!is.vector(slot(object, "xavini")))
        {
            return("slot >xavini< must be a vector!")
        }
        else if (!is.matrix(slot(object, "ini")))
        {
            return("slot >ini< must be a matrix!")
        }
        else if (!is.vector(slot(object, "Psi")))
        {
            return("slot >Psi< must be a vector!")
        }
        else  if (!is.matrix(slot(object, "lapla")))
        {
            return("slot >lapla< must be a matrix!")
        }


    }

 )
