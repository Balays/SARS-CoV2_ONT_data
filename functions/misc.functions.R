
dup <- function (x) x[duplicated(x)]

luniq <- function (x) length(unique(x))

## add prime 3s and prime 5s based on start, end and strand
add.primes <- function(df) {
    
  
  df$prime5[df$strand == '+'] <- df$start[df$strand == '+']
  df$prime3[df$strand == '+'] <- df$end[df$strand == '+'  ]
  ##
  df$prime5[df$strand == '-'] <- df$end[df$strand == '-']
  df$prime3[df$strand == '-'] <- df$start[df$strand == '-'  ]
  
  return(df)
  
}


