### A function that calculates the cross products of orthogonal polynomials
### up to a certain order.
#' @importFrom polynom polynomial
#' @export
Her = function(X_s, deg, type){
  ### Computing Hermite polynomials
  ### These are the terms in each variable up to the specified degree, without cross-products
  ### IF A FACTOR IS GIVEN AS THE SECOND ENTRY OF X_s, THEN IT COMPUTES THE CORRESPONDING DESIGN
  ### MATRIX BY REPEATING THE COMPUTATIONS THE NUMBER OF UNIQUE VALUES IN THE FACTOR.
  Her = c()
  if(!(dim(X_s)[2] %in% c(1,2))){print('Matrix X_s must have one or two columns. The second column
                                       must be a discrete variable')}
  else if(dim(X_s)[2] == 1){
    Her1 = polys(X_s[,1], which = type[1], normalized = TRUE, n = deg)
    Her = cbind(Her, Her1)
  }else{
    ### NOTE: TYPE ONLY VARIES FROM ONE ITERATION TO THE OTHER IF THE CONTINUOUS COVARIATE IN X_s
    ### IS DRAWN FROM DISTRIBUTIONS OF DIFFERENT SUPPORTS DEPENDING ON THE VALUE OF
    ### THE DISCRETE COVARIATE

    # COMPUTE MODEL MATRIX FROM FACTOR IN X_s
    factor = as.factor(X_s[,2])
    dummy = model.matrix(~ factor - 1)
  for(i in 1:length(unique(X_s[,2]))){
    Her1 = polys(X_s[,1], which = type[i], normalized = TRUE, n = deg)*dummy[, i]
    #for (k in 1:(deg+1)){
    #Remember the factor of 5 here
    #Her1[, k] = 5*EQL::hermite(X_s[,i], k-1, prob=FALSE)/sqrt((sqrt(pi)*2^(k-1)*factorial(k-1)))
    #}

    Her = cbind(Her, Her1)
  }
    temp = c()
    for(i in 1:(deg + 1)){
      temp = cbind(temp, Her[,i], Her[,i+deg+1])
    }
    Her = temp
  }

  # WE WONT USE THIS FOR NOW. ONLY CROSS PRODUCTS OF DISCRETE WITH CONTINUOUS VARIABLES FOR NOW
  if(0==1){
  ###Generating cross-products
  if(dim(X_s)[2]>1){
    Her2 = c()
    for(j in 1:dim(X_s)[1]){
      prov = Her[j, 1:(1+deg)]
      for (i in seq(from=1+deg+1,to=(dim(Her)[2]-deg), by = deg+1)){
        prov = prov%o%Her[j, i:(i+deg)]
      }
      ind = as.matrix((which(prov == prov,arr.ind = TRUE) -1))
      ind = ind[which(rowSums(ind)<=deg),]
      prov = prov[as.matrix(ind)+1]
      Her2=rbind(Her2,prov)
    }
    Her = Her2
  }
  }


  return(Her)
}


### Function to evaluate the orthogonal polynomial
###

polys = function(x, which, normalized = TRUE, n){
  if(which == "Laguerre"){
    leg4coef <- orthopolynom::laguerre.polynomials(n=n, normalized=normalized)
    leg4 <- as.matrix(as.data.frame(orthopolynom::polynomial.values(polynomials=leg4coef,
                                                                    x=x)))
  }
  if(which == "Hermite"){
    leg4coef <- orthopolynom::hermite.h.polynomials(n=n, normalized=normalized)
    leg4 <- as.matrix(as.data.frame(orthopolynom::polynomial.values(polynomials=leg4coef,
                                                                    x=x)))
  }
  if(which == "Legendre"){
    leg4coef <- orthopolynom::legendre.polynomials(n=n, normalized=normalized)
    leg4 <- as.matrix(as.data.frame(orthopolynom::polynomial.values(polynomials=leg4coef,
                                                                    x=x)))
  }

  if(which == "orth"){
    leg4 <- cbind(rep(1,length(x)), stats::poly(x,degree=n))
  }

  if(which == "Chebyshev"){
    leg4coef <- orthopolynom::chebyshev.s.polynomials(n=n, normalized=normalized)
    leg4 <- as.matrix(as.data.frame(orthopolynom::polynomial.values(polynomials=leg4coef,
                                                                    x=x)))
  }

  colnames(leg4) = 1:(n+1)
  return(leg4)
}

