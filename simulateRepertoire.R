# function to simulate activity spaces

simulate_aspace = function(
  part_codes, Locations, dist, gps_error, participants,
  aspaceSize, timeHome, where, propLoc, freqDurn,
  aspaceSizeModel, propLocModel, whereModel, fdModel, dist_bins)
{
  # with all model components now backfit location_types will be granular
  # and each model component picks the appropriate (potentially compound)
  # value, finally add FI for propLoc use
  location_types = c(as.character(unique(Locations$location_type)), "FI")
  location_types = location_types[-(which(location_types %in% c('ALOJAMIENTO','CEMENTERIOS')))]

  # sample activity space sizes
  if(is.null(dim(participants))){
    num_participants = 1
  } else {
    num_participants = nrow(participants)
  }
  if(aspaceSizeModel == 'Negative binomial')
    participants$aspace_size = rnbinom(num_participants, aspaceSize$nb$par[1], aspaceSize$nb$par[2])
  if(aspaceSizeModel == 'Poisson')
    participants$aspace_size = rpois(num_participants, aspaceSize$pois$par)
  if(aspaceSizeModel == 'Geometric')
    participants$aspace_size = rgeom(num_participants, aspaceSize$geom$par)

  # pull the correct propLoc parameters in the correct order of location_types names
  pl = propLoc[[propLocModel]]
  plParams = sapply(location_types, function(lt) {
    # find the location_type in propLoc
    par_index = grep(lt, names(pl$par))
    par = pl$par[par_index]
  })
  names(plParams) = location_types
  plParams = unlist(plParams)

  # sample location type distribution using multinomial
  type_counts = sapply(participants$aspace_size, function(size) {
    rmultinom(1, size, plParams)
  })
  type_counts = t(type_counts)
  
  # combine participants with their type counts
  type_counts = data.frame(type_counts)
  names(type_counts) = location_types
  participants = cbind(participants, type_counts)

  # indices to help select into the distance matrix by location_type
  indices_by_type = lapply(location_types, function(lt) {
    loc_ix = which(Locations$location_type==lt)
    return(loc_ix)
  })
  names(indices_by_type) = location_types
  
  # prep coefficients for all freq/durn parameters ignoring LL, k, and AICc
  fd = freqDurn[[fdModel]]
  fd = fd[1:(length(fd)-3)]
  
  # create appropriately named default values for each type
  default_values = c(1, 1, 0, 1, 1, 1, 0, 1, 0)
  names_values = c("freq.mean.a", "freq.mean.b", "freq.mean.c", "freq.std",
                    "durn.mean.a", "durn.mean.b", "durn.mean.c", "durn.std",
                    "corr")
  names(default_values) = names_values
  
  # pull best values from each type replacing neccesary default values
  # note we iterate to length(fd)-3 to avoid processing "LL", "k", and "AICc"
  fd_best = lapply(1:length(fd), function(i) {
    # start with appropriately named default values for this type
    type = fd[[i]]
    type_name = names(fd[i])
    type_values = default_values
    names(type_values) = paste0(names(type_values), ".", type_name)
    
    # and replace type_values with any params that actually exist in best
    best_values = type$best$par
    best_names = names(best_values)
    type_values[best_names] = best_values
    
    # but for use downstream move back to default names to avoid problems
    # with compound types
    names(type_values) = names_values
    return(type_values)
  })
  names(fd_best) = names(fd)
  

  # sample actual locations along with frequency and duration
  activity_spaces = lapply(1:num_participants, function(i)
  {
    # generate activity space for individual
    home_code = participants$home_code[i]
    aspace_ix = lapply(location_types[-which(location_types == 'FI')], function(lt) {
      # find the number of locations of this type to sample
      num_locs = participants[i, lt]

      # if the number of locations of this type is either 0 or the type is FI (fuera Iquitos, or outside Iquitos)
      # don't sample any locations
      if(num_locs==0) { return(numeric(0)) }

      # find distances to all locations of this type
      if(!is.null(dim(dist))){
        d = dist[home_code, indices_by_type[[lt]]]
      } else {
        d = dist[indices_by_type[[lt]]]
      }
      
      # find appropriate parameter with the where estimates
      where_model = where[[whereModel]]
      param_ix = grep(lt, names(where_model))
      param = where_model[[param_ix]]$par
      if(length(param) == 2){
        p = dexp(d ^ param[2], param[1])
      } else {
        p = dexp(d, param)
      }
      p = normalize(p)

      # sample location indices
      ix = sample(indices_by_type[[lt]], num_locs, replace=FALSE, p)
    })
    names(aspace_ix) = location_types[-which(location_types == 'FI')]
  
    # find distance and code for each activity space location
    aspace_dists = lapply(aspace_ix, function(ix)
      if(!is.null(dim(dist))){
        unname(dist[home_code, ix])
      } else {
        unname(dist[ix])
      })
    aspace_codes = lapply(aspace_ix, function(ix) as.character(Locations$location_code[ix]))
    aspace_codes = unlist(aspace_codes, use.names=FALSE)
    
    # sample frequencies to activity space locations
    aspace_freq_durn = lapply(location_types[-which(location_types == 'FI')], function(lt) {
      # number of locations of this type
      num_locs = length(aspace_dists[[lt]])
      if(num_locs==0) { return(numeric(0)) }
      
      # find corresponding best set for this type and pull params
      best = fd_best[[grep(lt, names(fd_best))]]
      freq.mean.a = best["freq.mean.a"]
      freq.mean.b = best["freq.mean.b"]
      freq.mean.c = best["freq.mean.c"]
      freq.std = best["freq.std"]
      durn.mean.a = best["durn.mean.a"]
      durn.mean.b = best["durn.mean.b"]
      durn.mean.c = best["durn.mean.c"]
      durn.std = best["durn.std"]
      
      # calculate log normal parameters for freq
      meanlog = .5 * freq.mean.a + freq.mean.a * freq.mean.b * exp(freq.mean.c * aspace_dists[[lt]]) /
        (exp(freq.mean.c * aspace_dists[[lt]]) + 1)
      sdlog = freq.std
      
      # sample frequencies
      freq = rlnorm(num_locs, meanlog, sdlog)

      # calculate log normal parameters for durn
      meanlog = .5 * durn.mean.a + durn.mean.a * durn.mean.b * exp(durn.mean.c * aspace_dists[[lt]]) /
        (exp(durn.mean.c * aspace_dists[[lt]]) + 1)
      sdlog = durn.std
      
      # sample durations
      durn = rlnorm(num_locs, meanlog, sdlog)
      
      return(list(freq=freq, durn=durn))
    })

    # add frequencies and durations for out of town locations
    numFI = participants[i, 'FI']
    if(numFI > 0){
      parFI = c('freq.mean.a.FI' = 0, 'freq.std.FI' = 1, 'durn.mean.a.FI' = 0, 'durn.std.FI' = 1, 'corr' = 0)
      parFI[names(freqDurn$FI$best$par)] = freqDurn$FI$best$par
      fdFI = exp(rmvnorm(
        numFI,
        parFI[c('freq.mean.a.FI', 'durn.mean.a.FI')],
        matrix(c(parFI['freq.std.FI']^2,
          prod(parFI[c('freq.std.FI', 'durn.std.FI', 'corr')]),
          prod(parFI[c('freq.std.FI', 'durn.std.FI', 'corr')]),
          parFI['durn.std.FI']^2), 2, 2)))

      aspace_freq_durn$FI = list(freq = fdFI[, 1], durn = fdFI[, 2])
      aspace_ix$FI = 1 : numFI
    }
    if(numFI == 0){
      aspace_freq_durn$FI = numeric()
      aspace_ix$FI = numeric()
    }

    # split into separate lists
    aspace_freqs = lapply(aspace_freq_durn, function(type) {
      if(length(type)==0) return(numeric(0))
      return(type$freq)
    })
    aspace_durns = lapply(aspace_freq_durn, function(type) {
      if(length(type)==0) return(numeric(0))
      return(type$durn)
    })
    
    # sample frequency and duration at home
    home_freq_durn =
      exp(rmvnorm(1, coef(timeHome$D)[c("freq.mean", "durn.mean")],
              cbind(c(coef(timeHome$D)["freq.std"]^2,
                      prod(coef(timeHome$D)[c("freq.std", "durn.std", "corr")])),
                    c(prod(coef(timeHome$D)[c("freq.std", "durn.std", "corr")]),
                      coef(timeHome$D)["durn.std"]^2))))
    home_freq = home_freq_durn[1, "freq.mean"]
    home_durn = home_freq_durn[1, "durn.mean"]
    
    # complete activity space from home and the sampled aspaces
    freq = c(home_freq, unlist(aspace_freqs, use.names=FALSE))    
    dur = c(home_durn, unlist(aspace_durns, use.names=FALSE))
    
    # step 6, generate Q matrix for each individual and final time allocation
    num_locs = length(freq)
    Q = matrix(0, num_locs, num_locs)
    for(ii in 1:num_locs) {
      for(jj in 1:num_locs) {
        if(ii == jj) {
          Q[ii, jj] = - 1 / dur[ii]
        }
        if(ii != jj) {
          Q[ii, jj] = (freq[jj] / sum(freq[-ii])) * (1 / dur[ii])
        }
      }
    }
    time.alloc = eigen(t(Q))$vectors[, num_locs]
    time.alloc = abs(time.alloc) / sum(abs(time.alloc))    
  
    # vector to index along locs, coords, dists, etc. with types for downstream
    aspace_types = sapply(location_types, function(type) {
      count = length(aspace_ix[[type]])
      if(count > 0 ) aspace(type, count)
      else c(as.character())
    })
    aspace_types = unlist(aspace_types, use.names=FALSE)
    aspace_types = c('home', aspace_types)
    
    list(types=aspace_types,
         time=time.alloc,
         codes=c(home_code,aspace_codes,rep("FI",sum(aspace_types=="FI"))))
  })
}
