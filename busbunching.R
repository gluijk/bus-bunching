# Bus bunching simulator
# www.overfitting.net
# https://www.overfitting.net/

library(ggplot2)
library(Cairo)


bus_bunching <- function(
        # headway: elapsed time between two consecutive bus arrivals at the same stop,
        # regardless of which bus arrives
        route_type = c("circular", "linear"),  # type of route
        n_buses = 6,              # number of buses in the simulation
        n_stops = 20,             # number of stops on the route
        route_length = 10000,     # total route length in m
        sim_time = 3*3600,        # total simulation time in s
        demand_rate = 0.05,       # passenger arrival rate per stop (passengers per second)
        # NOTE: at first passenger queues at stops are small or empty
        boarding_rate = 1.0,      # boarding rate (passengers per second)
        fixed_dwell = 5,          # fixed dwell time at each stop in s
        base_speed = 12,          # mean speed between stops (m/s)
        travel_sd = 2,            # standard deviation of travel time between stops in s
        control = c("none", "schedule", "headway"),  # optional holding control
        # Holding control policies applied at stops to mitigate bus bunching
        # Both "headway" and "schedule" operate by adding extra dwell time when a bus arrives "too early"
        # but differ in what timing reference they try to regulate
        # "none": fully unstable system, bunching emerges naturally
        # "headway": stabilizes service spacing locally, minimizes bunching but can cause drifting
        # "schedule": enforces nominal spacing, but because it still references observed headways it
        #             behaves like headway control with fixed target rather than true timetable adherence
        schedule_headway = NULL,  # target headway for schedule holding in s
        hold_limit = 60,          # maximum holding time when control is applied in s
        warmup = 0,               # warmup time to discard initial transient (0 = none) in s
        initial_headway = NULL,   # initial time delay between bus departures in s
        seed = NULL               # random seed for reproducibility

) {
    if (control == "schedule" && is.null(schedule_headway)) {
        message("ERROR: schedule_headway must be provided when control = 'schedule'")
        return(NULL)
    }
    
    if (!is.null(seed)) set.seed(seed)
    # match.arg() automatically picks the first element of the vector as the default:
    route_type <- match.arg(route_type)  # 'circular' by default
    control <- match.arg(control)  # 'none' by default
    
    # Uniform spacing assumption: only place where spatial structure enters the simulation
    seg_count <- n_stops  # number of travel segments
    seg_length <- route_length / seg_count  # distance between consecutive stops
    baseline_travel_time <- seg_length / base_speed  # mean travel time between stops 
    
    # Set default initial headway if not provided
    if (is.null(initial_headway)) {
        if (route_type == "circular") {
            loop_time_mean <- (baseline_travel_time + fixed_dwell) * seg_count
            initial_headway <- loop_time_mean / n_buses
        } else {
            initial_headway <- 0
        }
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
        new_arrivals <- rpois(1, demand_rate * delta_t)
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
                headway_obs <- t_now - last_arrival_time
                
                desired_headway = if (control == "schedule") schedule_headway else initial_headway
                hold_needed <- desired_headway - headway_obs
                if (hold_needed > 0) dwell_time <- dwell_time + min(hold_needed, hold_limit)
            }
        }
        
        # Travel implementation
        travel_time <- rnorm(1, baseline_travel_time, travel_sd)
        if (!is.finite(travel_time) || travel_time < 0.2)  travel_time <- baseline_travel_time
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
            mean_headway = mean_headway,
            cv_headway = cv_headway,
            mean_waiting_est = mean_waiting_est,
            percent_bunched = percent_bunched,
            n_events = nrow(events_df)  # total number of bus arrivals recorded after warm up filtering
            # For stable statistics it is desirable that n_events >> n_stops × n_buses
        )
    }
    
    # Output data
    list(
        events = events_df,  # dataframe
        headways = headways_df,  # datafame
        metrics = metrics  # simulation performance parameters
    )
}



################################################
# SIMULATION

# Minimal run (no control):
route_length=12000
demand_rate=0.02
boarding_rate=0.5  # 1.0
initial_headway = 60*8
schedule_headway = initial_headway*2
warmup=0

control="none"
control="headway"

for (control in c("none", "headway", "schedule")) {
    res <- bus_bunching(
        n_buses = 4, n_stops = 20,
        route_length = route_length,
        sim_time = 4 * 3600,
        demand_rate = demand_rate,
        boarding_rate = boarding_rate,  # NOTE: at first passenger queues at stops are small or empty
        base_speed = 10,
        travel_sd = 3,
        control = control,
        schedule_headway = schedule_headway, hold_limit = 60,
        route_type = "circular",
        initial_headway = initial_headway,
        warmup = warmup,
        seed = 1000
    )
    
    # res$metrics
    # head(res$events)
    
    
    # PNG output with antialiasing
    # Plot bus trajectories (arrival time vs stop index)
    CairoPNG(paste0("busbunching_", control, ".png"), width = 1920/2, height = 400)
        ev <- res$events
        p=ggplot(ev, aes(x = arrival/3600, y = stop, group = bus, color = factor(bus))) +
            geom_line(linewidth = 0.8) +
            geom_point(size = 1.2) +
            labs(
                x = "Time (hours)",
                y = "Bus stop",
                color = "Bus",
                title = paste0("Bus trajectories with control='", control, "' ",
                               ifelse(control == 'schedule', paste0("(schedule_headway=", round(schedule_headway/60,1), "min) "), ""),
                               "(initial_headway=", round(initial_headway/60,1),"min, route_length=", route_length,
                               "m, demand_rate=", demand_rate, "pax/s, boarding_rate=",
                               boarding_rate, "pax/s)")
            ) +
            theme_minimal(base_size = 22) +
            theme(
                legend.title = element_text(size = 20),
                legend.text  = element_text(size = 20),
                axis.title   = element_text(size = 20),
                axis.text    = element_text(size = 20),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_line(linewidth = 0.6),
                axis.line = element_line(linewidth = 1),
                legend.key.size = unit(1.5, "lines"),
                plot.title = element_text(size = 14)  # reduce title size here
            )
        
    print(p)
    dev.off()
}


################################################
# PLOTTING

# Plotting headways over time
hw <- res$headways
png("busbunchingnone_headways.png", width=960, height=400)
ggplot(hw, aes(x = arrival/3600, y = headway/60)) +  # global aesthetics
    geom_point(alpha = 0.6) +
    geom_smooth(se = FALSE) +  # linear regresion (~headway average)
    labs(x = "Time (hours)", y = "Headway (minutes)",
         title = "Headways at stops over time") +
    scale_y_continuous(limits = c(0, 30)) +  # set y-axis maximum to 30
    theme_minimal()
dev.off()

# Headway distribution
png("busbunchingnone_hist.png", width=960, height=280)
hist(hw$headway/60, breaks=150, xlab='Headway (minutes)', xlim=c(0,30), main='Headway distribution')
abline(v=c(mean(hw$headway/60), median(hw$headway/60)), lty='dotted', col='red')
dev.off()

# Plot bus trajectories (arrival time vs stop index)
png("busbunchingpepe.png", width=1920, height=1080)
    ev <- res$events
    ggplot(ev, aes(x = arrival/3600, y = stop, group = bus, color = factor(bus))) +
        geom_line() + geom_point(size = 0.6) +
        labs(x = "Time (hours)", y = "Stop index", color = "Bus") +
        theme_minimal()
dev.off()


