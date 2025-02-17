

mixest2 <- function(y,x,mods=NULL,ftype=NULL,V=NULL,W=NULL,atype=NULL,Tvar=NULL, max_vars = NULL)
  {
    if (is.null(atype)) { atype <- 0 }
    
    if (atype==0)
      {
        out <- .mixest2a(y=y,x=x,mods=mods,ftype=ftype,V=V,W=W,Tvar=Tvar,max_vars = max_vars)
      }
    else
      {
        out <- .mixest2b(y=y,x=x,mods=mods,ftype=ftype,V=V,W=W,Tvar=Tvar,max_vars = max_vars)
      }

    return(out)
  }
  
  