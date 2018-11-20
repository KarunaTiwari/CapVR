#' Vehicle routing using local search
#'
#' A local search algorithm is implemented to determine the optimal routes.
#' The result can be fed into the Genetic algorithm based routing function for further optimisation.
#' @param loc_data A data frame with the nodal points and the corresponding number of people to be dropped/picked up from the nodes.
#' @param min_cap The minimum number of people that should be allotted to each vehicle.
#' @param dist_mat The distance matrix to be used. Dimensions should be Px(P+1), where P is the number of nodes.
#' @param time_mat The time matrix to be used. Dimensions should be Px(P+1), where P is the number of nodes.
#' @param loc_code The column number/name in dist_mat which contains the nodal points.
#' @param init_node The name of the starting point/node.
#' @return A list containing tables for each cab allocation along with distance and time computed.
#' The list also contains a vectorised solution which can be fed into the GA based routing algorithm.
#' @export

LS_Route<-function(loc_data,loc_code,min_cap,dist_mat,time_mat,init_node){

  require(data.table)

  dt<-as.data.table(loc_data)
  colnames(dt)<-c("PickupCode","Employees")
  dist_mat<-as.data.table(dist_mat)
  time_mat<-as.data.table(time_mat)

  if(class(loc_code) == "numeric"){
    pick<-loc_code
  }
  else if(class(loc_code) == "character"){
    pick<-which(colnames(dist_mat) == loc_code)
  }
  ori<-dist_mat[[pick]]
  dist_mat<-dist_mat[,-pick,with = F]
  time_mat<-time_mat[,-pick,with = F]

  rem_pick<-dt$PickupCode[which(dt$Employees <= 0 & dt$PickupCode != init_node)]
  rem_mat<-which(colnames(dist_mat) %in% rem_pick)
  dist_mat<-dist_mat[-(rem_mat),-(rem_mat),with = F]
  time_mat<-time_mat[-(rem_mat),-(rem_mat),with = F]
  dt<-dt[Employees > 0]

  res_list<-vector(mode = "list")
  r<-1

  while(sum(dt$Employees) > 0){
    init<-which(colnames(dist_mat) == init_node)
    cab_cnt<-0
    w<-1
    stop_dist<-emp_add<-stop_time<-numeric()
    while(cab_cnt < min_cap & sum(dt$Employees) > 0){
      visit<-as.numeric(dist_mat[init,])[-init]
      names(visit)<-colnames(dist_mat)[-init]
      visit_t<-as.numeric(time_mat[init,])[-init]
      names(visit_t)<-colnames(time_mat)[-init]
      if(w == 1){
        stop<-names(visit)[which.min(visit)]
        stop_dist<-c(stop_dist,visit[stop])
        stop_time<-c(stop_time,visit_t[stop])
      }
      else{
        stop<-names(visit[-which(names(visit) == init_node)])[which.min(visit[-which(names(visit) == init_node)])]
        stop_dist<-c(stop_dist,(visit[stop] + stop_dist[length(stop_dist)]))
        stop_time<-c(stop_time,(visit_t[stop] + stop_time[length(stop_time)]))
      }
      emp_cnt<-dt[PickupCode == stop]$Employees[1]
      if(is.na(emp_cnt) == T){
        emp_cnt<-0
      }
      if((min_cap - cab_cnt) <= emp_cnt){
        if(w == 1){
          emp_add<-c(emp_add,(min_cap - cab_cnt))
          dt[PickupCode == stop,Employees := Employees - (min_cap - cab_cnt)]
          cab_cnt<-min_cap
        }
        else {
          emp_add<-c(emp_add,(min_cap - cab_cnt))
          dt[PickupCode == stop,Employees := Employees - (min_cap - cab_cnt)]
          cab_cnt<-min_cap
          dist_mat<-dist_mat[-init,-init,with = F]
          time_mat<-time_mat[-init,-init,with = F]
          dt<-dt[Employees > 0]
        }
      }
      else {
        if(w == 1){
          cab_cnt<-cab_cnt + emp_cnt
          emp_add<-c(emp_add,emp_cnt)
          dt[PickupCode == stop,Employees := Employees - emp_cnt]
          dt<-dt[Employees > 0]
          init<-which(colnames(dist_mat) == stop)
        }
        else {
          cab_cnt<-cab_cnt + emp_cnt
          emp_add<-c(emp_add,emp_cnt)
          dt[PickupCode == stop,Employees := Employees - emp_cnt]
          dist_mat<-dist_mat[-init,-init,with = F]
          time_mat<-time_mat[-init,-init,with = F]
          dt<-dt[Employees > 0]
          init<-which(colnames(dist_mat) == stop)
        }
      }
      w<-w+1
    }
    res_int<-data.table("Node" = names(stop_dist),"Cum_Distance" = stop_dist,"Cum_Time" = stop_time,"Employees" = emp_add)
    res_list[[r]]<-res_int
    r<-r+1
  }

  names(res_list)<-paste("Cab",1:length(res_list),sep = "-")
  res_list_2<-lapply(X = res_list,FUN = function(x){
    zer<-which(x$Employees == 0)
    if(is.na(zer) == F && length(zer) >= 1){
      res<-x[-zer,]
    }
    else {
      res<-x
    }
    return(res)
  })

  res_3<-unlist(lapply(X = res_list_2,FUN = function(x){
    tour<-character()
    for(i in 1:nrow(x)){
      tour<-c(tour,rep(x$Node[i],times = x$Employees[i]))
    }
    return(tour)
  }))
  names(res_3)<-NULL

  ress<-list("Tables" = res_list_2,"Vector" = res_3)

  return(ress)

}

#' Vehicle routing using Genetic algorithm.
#'
#' GA is implemented to determine the best routes so as to minimise the total distance covered by all vehicles.
#' The algorithm converges much faster and presents better results when the vectorised output from the LS_Route function
#' is provided as a suggested solution.
#' @param loc_data A data frame with the nodal points and the corresponding number of people to be dropped/picked up from the nodes.
#' @param dist_mat The distance matrix to be used. Dimensions should be Px(P+1), where P is the number of nodes.
#' @param min_cap The minimum number of people that should be allotted to each vehicle.
#' @param init_node The name of the starting point/node.
#' @param sug_sol Optional suggested solution to initiate the GA. The vectorised output from the LS_Route function should ideally be used.
#' @param mxiter Maximum number of iterations for which the GA should be run. 50000 by default.
#' @param mut_rate The rate of mutation to be used. 0.2 by default.
#' @param crs_rate The rate of crossover to be used. 0.8 by default.
#' @param converge The maximum number of iterations to be run before stopping the GA if there is no improvement in the solution. 5000 by default.
#' @param re_eval Should the GA solution be re-evaluated for further optimisation.
#' @return A list with the optimised vehicle routes and the total distance.
#' @export

GA_Route<-function(loc_data,dist_mat,min_cap,init_node,sug_sol = NULL,mxiter = 50000,mut_rate = 0.2,crs_rate = 0.8,converge = 5000,re_eval = T){

  require(data.table)
  require(GA)

  dt_pick<-as.data.table(loc_data)
  colnames(dt_pick)<-c("PickupCode","Employees")
  D<-as.data.frame(dist_mat)
  rownames(D)<-D[,1]
  D<-D[,-1]

  P<-character()
  for(i in 1:nrow(dt_pick)){
    P<-c(P,rep(dt_pick$PickupCode[i],times = dt_pick$Employees[i]))
  }

  dist_opt<-function(x){
    ls<-split(x, ceiling(seq_along(x)/min_cap))
    ls<-lapply(X = ls,FUN = function(x){
      P[x]
    })
    lsd<-lapply(X = ls,FUN = function(x){
      n<-init_node
      d<-numeric()
      for(i in 1:length(x)){
        d[i]<-D[n,x[i]]
        n<-x[i]
      }
      return(sum(d))
    })
    dist<-unlist(lsd)
    return(sum(dist))
  }

  fitness<-function(x){
    -dist_opt(x)
  }

  if(is.null(sug_sol) == T){
    mod_ga<-ga(type = "permutation",fitness = fitness,lower = 1,upper = length(P),maxiter = mxiter,
               pmutation = mut_rate,pcrossover = crs_rate,run = converge,monitor = plot)
  }

  else {
    sol<-numeric()
    for(i in 1:length(sug_sol)){
      pos<-which(P == sug_sol[i])
      sol[i]<-pos[which(pos %in% sol == F)[1]]
    }

    mod_ga<-ga(type = "permutation",fitness = fitness,lower = 1,upper = length(P),maxiter = mxiter,monitor = T,
               pmutation = mut_rate,pcrossover = crs_rate,run = converge,suggestions = matrix(data = sol,ncol = length(sol),nrow = 1))
  }

  dist<-function(x){
    ls<-split(x, ceiling(seq_along(x)/min_cap))
    names(ls)<-paste("Cab",1:length(ls),sep = "-")
    lsd<-lapply(X = ls,FUN = function(x){
      n<-init_node
      d<-numeric()
      for(i in 1:length(x)){
        d[i]<-D[n,x[i]]
        n<-x[i]
      }
      return(sum(d))
    })
    tot_dist<-unlist(lsd)
    return(list("Solution" = ls,"Distance" = sum(tot_dist)))
  }

  ga_sol<-dist(P[as.numeric(mod_ga@solution[1,])])

  if(re_eval == F){
    return(ga_sol)
  }

  else {
    multi_ga<-function(l){
      route<-l
      tryCatch(expr = {
        opt_cab<-function(x){
          n<-init_node
          cab_dist<-numeric()
          x<-route[x]
          for(i in 1:length(x)){
            cab_dist[i]<-D[n,x[i]]
            n<-x[i]
          }
          return(sum(cab_dist))
        }
        fitness_2<-function(x){
          -opt_cab(x)
        }
        mod_cab<-ga(type = "permutation",fitness = fitness_2,lower = 1,upper = length(route),
                    pcrossover = crs_rate,pmutation = mut_rate,maxiter = 1000,
                    suggestions = matrix(data = 1:length(route),nrow = 1,ncol = length(route)))
        sol<-route[as.numeric(mod_cab@solution[1,])]

        return(sol)
      },error = function(e){
        return(l)
      })
    }

    res<-lapply(X = ga_sol$Solution,FUN = multi_ga)
    res_vec<-(unlist(res))
    names(res_vec)<-NULL
    fin_sol<-dist(res_vec)

    return(fin_sol)

  }
}
