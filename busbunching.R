# Basic bus bunching discrete event simulation
# www.overfitting.net
# https://www.overfitting.net/2025/12/simulador-de-bus-bunching-con-r.html

library(ggplot2)
library(Cairo)


# Discrete Event Simulation (DES) for Bus bunching
# Modeling approach in which the system state changes only at a discrete set of time points called events
# DES models are: order sensitive, causality driven, state based
bus_bunching <- function(

    # Route definition:
    route_type = c("circular", "linear"),  # type of route
    route_length = 10000,     # total route length in m
    n_stops = 20,             # number of stops on the route
    n_buses = 6,              # number of buses in the simulation
    fixed_dwell = 5,          # fixed dwell time at each stop in s
    initial_headway = NULL,   # initial time delay between bus departures in s
    # headway: elapsed time between two consecutive bus arrivals at the same stop,
    # regardless of which bus arrives
    base_speed = 12,          # mean speed between stops in m/s
    travel_sd_fraction = 0.2, # sd of travel time between stops as a fraction of baseline_travel_time
    
    # Demand definition:
    # NOTE: at first passenger queues at stops are small or empty -> first bus rides fast
    demand_rate = 0.05,       # passenger arrival rate per stop in passengers/second
    boarding_rate = 1.0,      # boarding rate in passengers/second
    
    # Holding control:
    control = c("none", "schedule", "headway"),
    # Policies applied at stops to mitigate bus bunching. Both "headway" and "schedule"
    # operate by adding extra dwell time when a bus arrives "too early"
    # but differ in what timing reference they try to regulate
    # "none": fully unstable system, bunching emerges naturally
    # "headway": stabilizes service spacing locally, minimizes bunching but can cause drifting
    # "schedule": enforces nominal spacing, but because it still references observed headways it
    #             behaves like headway control with fixed target rather than true timetable adherence
    schedule_headway = NULL,  # target headway for schedule holding in s
    hold_limit = 60,          # maximum holding time when control is applied in s
    
    # Simulation parameters:
    warmup = 0,               # warmup time to discard initial transient (0=none) in s
    sim_time = 3*3600,        # total simulation time in s
    seed = 123                # random seed for reproducibility
    
) {
    # Optional parameters: match.arg() automatically picks the first element of the vector as default
    route_type <- match.arg(route_type)  # 'circular' by default
    control <- match.arg(control)  # 'none' by default
    if (control == "schedule" && is.null(schedule_headway)) {
        message("ERROR: schedule_headway must be provided when control = 'schedule'")
        return(NULL)
    }
    if (!is.null(seed)) set.seed(seed)
    
    # Uniform spacing assumption: only place where spatial structure enters the simulation
    seg_count <- n_stops  # number of travel segments
    seg_length <- route_length / seg_count  # distance between consecutive stops
    baseline_travel_time <- seg_length / base_speed  # mean travel time between stops
    min_travel_time <- 0.01 * baseline_travel_time  # min travel_time to preserve event ordering (1% of baseline)
    travel_time_sd <- travel_sd_fraction * baseline_travel_time  # convert fractional SD to s
    
    # Set default initial headway if not provided
    if (is.null(initial_headway)) {
        loop_time_mean <- (baseline_travel_time + fixed_dwell) * seg_count
        initial_headway <- loop_time_mean / n_buses
        # This would be the exact initial_headway with no passengers (demand_rate=0) to
        # make the last bus depart when the first bus is about to finish the route
        # In a real application with passengers the ideal initial_headway should be set higher
    }
    
    # Initialize buses
    buses <- data.frame(
        bus = seq_len(n_buses),
        next_stop = rep(1, n_buses),
        next_arrival = seq(0, by = initial_headway, length.out = n_buses),
        stringsAsFactors = FALSE
    )
    
    last_visit_time <- rep(0, seg_count)
    waiting_passengers <- rep(0L, seg_count)
    events <- vector("list", 0)
    
    record_event <- function(bus_id, stop_id, t_arrival, waiting, boarded, dwell, travel_time) {
        events[[length(events) + 1L]] <<- list(
            # Unlike <- (regular assignment, which assigns in the current environment)
            # <<- assigns to a variable in the parent environment or higher
            bus = bus_id,
            stop = stop_id,
            arrival = t_arrival,
            waiting = waiting,
            boarded = boarded,
            dwell = dwell,
            travel_time = travel_time
        )
    }
    
    while (TRUE) {
        next_idx <- which.min(buses$next_arrival)
        t_now <- buses$next_arrival[next_idx]
        if (t_now > sim_time) break
        
        bus_id <- buses$bus[next_idx]
        stop_id <- buses$next_stop[next_idx]
        
        delta_t <- t_now - last_visit_time[stop_id]
        if (delta_t < 0) delta_t <- 0
        new_arrivals <- rpois(1, demand_rate * delta_t)  # Poisson arrivals
        waiting_passengers[stop_id] <- waiting_passengers[stop_id] + new_arrivals
        last_visit_time[stop_id] <- t_now
        
        boarded <- waiting_passengers[stop_id]
        waiting_passengers[stop_id] <- 0L
        dwell_time <- fixed_dwell + boarded / boarding_rate
        
        # Holding control
        if (control %in% c("schedule", "headway")) {
            prev_idx <- which(sapply(events, function(x) x$stop) == stop_id)
            if (length(prev_idx) > 0) {
                last_arrival_time <- events[[tail(prev_idx, 1)]]$arrival
                headway_obs <- t_now - last_arrival_time  # headway observed from last departed bus
                
                desired_headway = if (control == "schedule") schedule_headway else initial_headway
                hold_needed <- desired_headway - headway_obs
                if (hold_needed > 0) dwell_time <- dwell_time + min(hold_needed, hold_limit)
            }
        }
        
        # Travel implementation
        travel_time <- max(min_travel_time, rnorm(1, baseline_travel_time, travel_time_sd))
        record_event(bus_id, stop_id, t_now, new_arrivals, boarded, dwell_time, travel_time)
        next_stop <- stop_id + 1
        
        if (route_type == "circular") {
            if (next_stop > seg_count) next_stop <- 1
            buses$next_stop[next_idx] <- next_stop
            buses$next_arrival[next_idx] <- t_now + dwell_time + travel_time
        } else {  # linear
            if (next_stop > seg_count) {
                buses$next_arrival[next_idx] <- Inf
            } else {
                buses$next_stop[next_idx] <- next_stop
                buses$next_arrival[next_idx] <- t_now + dwell_time + travel_time
            }
        }
    }
    
    if (length(events) == 0) {
        events_df <- data.frame()
    } else {
        events_df <- do.call(rbind, lapply(events, function(e) {
            data.frame(
                bus = e$bus,
                stop = e$stop,
                arrival = e$arrival,
                waiting = e$waiting,
                boarded = e$boarded,
                dwell = e$dwell,
                travel_time = e$travel_time,
                stringsAsFactors = FALSE
            )
        }))
    }
    
    # Compute headways after warmup period
    if (nrow(events_df) > 0 && warmup > 0) {
        events_df <- subset(events_df, arrival > warmup)
    }
    
    
    # Headways dataframe calculation
    headways_list <- lapply(split(events_df, events_df$stop), function(dfstop) {
        if (nrow(dfstop) <= 1) return(NULL)
        times <- sort(dfstop$arrival)
        data.frame(
            stop = dfstop$stop[1],
            arrival = times[-1],
            headway = diff(times)
        )
    })
    headways_df <- if (length(headways_list) == 0) data.frame() else do.call(rbind, headways_list)
    
    
    # Metrics calculation
    if (nrow(headways_df) == 0) {
        metrics <- list(mean_headway = NA, cv_headway = NA, mean_waiting_est = NA, percent_bunched = NA)
    } else {
        hw <- subset(headways_df, arrival > warmup)
        
        mean_headway <- mean(hw$headway)  # average time between consecutive bus arrivals at a stop
        cv_headway <- sd(hw$headway) / mean_headway  # relative dispersion of headways around their mean
        # In the context of bunching this is the primary stability metric:
        #   0.0–0.2: very regular service, 0.2–0.4: mild irregularity
        #   0.4–0.6: noticeable bunching, >0.6: severe bunching
        
        mean_waiting_est <- mean(hw$headway) / 2  # estimated mean passenger waiting time at stops
        
        med_h <- median(hw$headway)  # median of headways
        # Fraction of headways abnormally long relative to typical service
        percent_bunched <- sum(hw$headway > 2 * med_h) / nrow(hw)
        #   <0.05: very stable, 0.05–0.15: mild bunching
        #   0.15–0.30: frequent bunching, >0.30: systemically unstable
        
        # cv_headway: measures dispersion globally
        # percent_bunched: measures operational failure frequency
        # Together they are the most informative pair
        
        metrics <- list(
            # For stable statistics it is desirable that n_events >> n_stops × n_buses
            n_events = nrow(events_df),  # total number of bus arrivals recorded after warm up filtering
            mean_headway = mean_headway,
            mean_waiting_est = mean_waiting_est,  # half of mean_headway
            cv_headway = cv_headway,
            percent_bunched = percent_bunched
        )
    }
    
    
    # Record input parameters actually used by the simulation
    params <- list(
        route_type = route_type,
        route_length = route_length,
        n_stops = n_stops,
        n_buses = n_buses,
        fixed_dwell = fixed_dwell,
        initial_headway = initial_headway,
        base_speed = base_speed,
        travel_sd_fraction = travel_sd_fraction,  # fractional SD input
        travel_time_sd = travel_time_sd,  # actual SD used in s
        demand_rate = demand_rate,
        boarding_rate = boarding_rate,
        control = control,
        schedule_headway = schedule_headway,
        hold_limit = hold_limit,
        warmup = warmup,
        sim_time = sim_time,
        seed = seed
    )
    
    
    # Output data
    list(
        events = events_df,  # dataframe
        headways = headways_df,  # datafame
        metrics = metrics,  # simulation performance parameters
        params = params  # input parameters actually used
    )
}


################################################

# 1. BASIC EXAMPLE: 3 BUSES LINE WITH DIFFERENT DEMAND AND CONTROL

demand_rate <- c(0.00, 0.01, 0.05, 0.05)
control <- c("none", "none", "none", "headway")
hold_limit <- c(NA, NA, NA, 60*2)

for (i in 1:4) {
    res <- bus_bunching(
        n_buses = 3, n_stops = 8,
        route_length = 15000,
        demand_rate = demand_rate[i],
        initial_headway = 60 * 15,
        base_speed = 5,
        sim_time = 4 * 3600,
        control = control[i],
        hold_limit = hold_limit[i]
    )
    
    res$metrics

    
    # Plot bus trajectories (arrival time vs stop index)
    ev <- res$events
    params <- res$params
    for (f in c(1,3)) {
        CairoPNG(paste0("busbunching", f, "_sim", i, ".png"), width = 512*f, height = 220*f)
            p=ggplot(ev, aes(x = arrival/3600, y = stop, group = bus, color = factor(bus))) +
                geom_line(linewidth = 0.8*f) +
                geom_point(size = 1.8*f) +
                labs(
                    x = "Time (hours)",
                    y = "Bus stop",
                    color = "Bus",
                    title = paste0("Bus trajectories with control='", params$control, "' ",
                       ifelse(params$control == 'schedule', paste0("[ schedule_headway=",
                       round(params$schedule_headway/60,1), "min ] "), ""),
                       "\n[ initial_headway=", round(params$initial_headway/60,1),"min / route_length=",
                       params$route_length, "m ]\n[ demand_rate=", params$demand_rate*60, "pax/min / boarding_rate=",
                       params$boarding_rate, "pax/s ]")
                ) +
                theme_minimal(base_size = 22) +
                theme(
                    legend.title = element_text(size = 12*f),
                    legend.text  = element_text(size = 12*f),
                    axis.title   = element_text(size = 12*f),
                    axis.text    = element_text(size = 12*f),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_line(linewidth = 0.6),
                    axis.line = element_line(linewidth = 1),
                    legend.key.size = unit(1.5, "lines"),
                    plot.title = element_text(size = 12*f)  # reduce title size here
                )
            
            print(p)
        dev.off()
    }
}


################################################

# 2. DENSE OPTIMIZED BUS LINE (TRAJECTORIES, HEADWAYS, HISTOGRAM)


control <- c("none", "schedule")
hheadway <- c(60 * 3, 60 * 4)

for (i in 1:2) {
    f=1
    res <- bus_bunching(
        n_buses = 10, n_stops = 20,
        route_length = 10000,
        demand_rate = 0.02 * f * 1.5,
        initial_headway = 60 * 3,
        schedule_headway = 60 * 4,  # on a regular basis we try to expand the headway to 4min
        travel_sd_fraction = 0.2,
        base_speed = 5,
        sim_time = 3600 * 5,
        control = control[i],
        hold_limit = 60,  # but with a max waiting time of 1min
        warmup = 3600 * 3.7 * 0  # testing the warmup parameter -> works fine
    )
    
    res$metrics
    
    
    # Plot bus trajectories (arrival time vs stop index)
    f=2
    ev <- res$events
    params <- res$params
    CairoPNG(paste0("busbunchingdense_", control[i], ".png"), width = 1920/2*f, height = 1080/4*f)
        p=ggplot(ev, aes(x = arrival/3600, y = stop, group = bus, color = factor(bus))) +
        geom_line(linewidth = 0.2*f) +
        geom_point(size = 0.6*f) +
        labs(
            x = "Time (hours)",
            y = "Bus stop",
            color = "Bus",
            title = paste0("Bus trajectories with control='", params$control, "' ",
               ifelse(params$control == 'schedule', paste0("[ schedule_headway=",
                                                           round(params$schedule_headway/60,1), "min ] "), ""),
               "\n[ initial_headway=", round(params$initial_headway/60,1),"min / route_length=",
               params$route_length, "m ]\n[ demand_rate=", params$demand_rate*60, "pax/min / boarding_rate=",
               params$boarding_rate, "pax/s ]")
        ) +
        theme_minimal(base_size = 22) +
        theme(
            legend.title = element_text(size = 8*f),
            legend.text  = element_text(size = 8*f),
            axis.title   = element_text(size = 8*f),
            axis.text    = element_text(size = 8*f),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_line(linewidth = 0.6),
            axis.line = element_line(linewidth = 1),
            legend.key.size = unit(1.5, "lines"),
            plot.title = element_text(size = 8*f)  # reduce title size here
        )
    
        print(p)
    dev.off()
    
    
    # Headways over time
    hw <- res$headways
    CairoPNG(paste0("busbunchingdense_headways_", control[i],".png"), width = 512, height = 400)
        plot(hw$arrival / 3600, hw$headway / 60,
             pch = 16, cex = 0.75, col = rgb(0, 0, 0, alpha = 0.3),
             xlab = "Time (hours)", ylab = "Headway (minutes)", main = "Headways at stops over time",
             ylim = c(0, 25),
             cex.lab = 1.2, cex.main = 1.3, cex.axis = 1.0
        )
        abline(h=hheadway[i] / 60, col = "red", lwd = 1)
    dev.off()
    
    
    # Headways distribution
    CairoPNG(paste0("busbunchingdense_hist_", control[i], ".png"), width=512, height=300)
        hist(hw$headway/60, breaks=500, xlab='Headway (minutes)',
             xlim=c(0, 25), main='Headway distribution')
        abline(v=c(mean(hw$headway/60), median(hw$headway/60)), lty='dotted', col='red')
    dev.off()
}


