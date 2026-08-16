# This script contains various similarity functions 
# used to perform Fuzzy Time Series (FTS) calculations.

library(DescTools)

Similarity <- function(x, y, E, Tn, U) {
  # Calculates the similarity between two intervals.
  # Parameters:
  # x, y: Intervals to compare
  # E: Embedding choice
  #   E=1: Length-based
  #   E=2: Lukasiewicz
  #   E=3: Gödel
  #   E=4: Goguen
  #   E=5: Fodor
  #   E=6: Exponential
  # Tn: T-norm 
  #   Tn=1: Arithmetic mean
  #   Tn=2: Minimum
  #   Tn=3: Geometric mean
  #   Tn=4: Maximum
  # U: Aggregation union
  #   U=0: Standard similarity (no union)
  #   U=1: Aggregated similarity
  
  if (E == 1) {
    Inc <- function(a, b) {
      # Calculates length-based embedding for intervals a and b
      a_low <- a[1]
      a_high <- a[2]
      b_low <- b[1]
      b_high <- b[2]
      
      # Check for intersection
      if (max(a_low, b_low) > min(a_high, b_high)) {
        return(0.0)
      }
      # Check if a is contained in b
      if (a_low >= b_low & a_high <= b_high) {
        return(1.0)
      }
      
      return(Overlap(a, b) / (a[2] - a[1]))
    }
    
  } else if (E == 2) {
    Inc <- function(a, b) {
      # Calculates Lukasiewicz embedding for intervals a and b
      a_low <- a[1]
      a_high <- a[2]
      b_low <- b[1]
      b_high <- b[2]
      
      # Check for intersection
      if (max(a_low, b_low) > min(a_high, b_high)) {
        return(0.0)
      }
      
      # General form
      term1 <- 1 - b_low + a_low
      term2 <- 1 - a_high + b_high
      return(min(max(term1, 0), max(term2, 0), 1))
    }
    
  } else if (E == 3) {
    Inc <- function(a, b) {
      # Calculates Gödel embedding for intervals a and b
      a_low <- a[1]
      a_high <- a[2]
      b_low <- b[1]
      b_high <- b[2]
      
      # Check for intersection
      if (max(a_low, b_low) > min(a_high, b_high)) {
        return(0.0)
      }
      
      # Check if a is contained in b
      if (a_low >= b_low & a_high <= b_high) {
        return(1.0)
      }
      
      # Case: b_low <= a_low <= b_high < a_high
      if (b_low <= a_low & a_low <= b_high & b_high < a_high) {
        return(b_high)
      }
      
      # Default case
      return(a_low)
    }
    
  } else if (E == 4) {
    Inc <- function(a, b) {
      # Calculates Goguen embedding for intervals a and b
      a_low <- a[1]
      a_high <- a[2]
      b_low <- b[1]
      b_high <- b[2]
      
      # Check for intersection
      if (max(a_low, b_low) > min(a_high, b_high)) {
        return(0.0)
      }
      
      # Check if a is contained in b
      if (a_low >= b_low & a_high <= b_high) {
        return(1.0)
      }
      
      # Case: a ∩ b ≠ ∅ and b_low = 0
      if (b_low == 0) {
        return(b_high / a_high)
      }
      
      # General case
      term1 <- b_high / a_high
      term2 <- a_low / b_low
      return(min(term1, term2))
    }
    
  } else if (E == 5) {
    Inc <- function(a, b) {
      # Calculates Fodor embedding for intervals a and b
      a_low <- a[1]
      a_high <- a[2]
      b_low <- b[1]
      b_high <- b[2]
      
      # Check for intersection
      if (max(a_low, b_low) > min(a_high, b_high)) {
        return(0.0)
      }
      
      # Check if a is contained in b
      if (a_low >= b_low & a_high <= b_high) {
        return(1.0)
      }
      
      # Case: b_low <= a_low <= b_high < a_high
      if (b_low <= a_low & a_low <= b_high & b_high < a_high) {
        return(max(1 - a_high, b_high))
      }
      
      # Case: a_low < b_low <= a_high <= b_high
      if (a_low < b_low & b_low <= a_high & a_high <= b_high) {
        return(max(1 - b_low, a_low))
      }
      
      # General case
      term1 <- max(1 - a_high, b_high)
      term2 <- max(1 - b_low, a_low)
      return(min(term1, term2))
    }
    
  } else if (E == 6) {
    Inc <- function(a, b) {
      # Calculates Exponential embedding for intervals a and b
      a_low <- a[1]
      a_high <- a[2]
      b_low <- b[1]
      b_high <- b[2]
      
      # Check for intersection
      if (max(a_low, b_low) > min(a_high, b_high)) {
        return(0.0)
      }
      
      # Check if a is contained in b
      if (a_low >= b_low & a_high <= b_high) {
        return(1.0)
      }
      
      # Case: b_low <= a_low <= b_high < a_high
      if (b_low <= a_low & a_low <= b_high & b_high < a_high) {
        return(b_high * exp(1 - a_high))
      }
      
      # Case: a_low < b_low <= a_high <= b_high
      if (a_low < b_low & b_low <= a_high & a_high <= b_high) {
        return(a_low * exp(1 - b_low))
      }
      
      # General case
      term1 <- b_high * exp(1 - a_high)
      term2 <- a_low * exp(1 - b_low)
      return(min(term1, term2))
    }
  }
  
  # T-norm 
  if (Tn == 1) {
    Tnorm <- function(a, b) { return((a + b) / 2) }
  } else if (Tn == 2) {
    Tnorm <- function(a, b) { return(min(a, b)) }
  } else if (Tn == 3) {
    Tnorm <- function(a, b) { return(sqrt(a * b)) }
  } else if (Tn == 4) {
    Tnorm <- function(a, b) { return(max(a, b)) }
  }
  
  # Aggregation union function
  aggregation <- function(a, b) {
    aub <- c(min(a[1], b[1]), max(a[2], b[2]))
    return(aub)
  }
  
  if (U != 0) {
    return(Tnorm(Inc(aggregation(x, y), x), Inc(aggregation(x, y), y)))
  } else {
    return(Tnorm(Inc(x, y), Inc(y, x)))
  }
}
