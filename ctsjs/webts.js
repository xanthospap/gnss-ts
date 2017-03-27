// ---------------------------------------------------------------------------
//  function deepCopy
//
//  Return a (deep) copy of the input array.
//  --------------------------------------------------------------------------
function deepCopy(arr)
{
    var out = [];
    for (var i = 0, len = arr.length; i < len; i++) {
        var item = arr[i];
        var obj  = {};
        for (var k in item) { obj[k] = item[k]; }
        out.push(obj);
    }
    return out;
}

// ---------------------------------------------------------------------------
//  function mjd_to_gps
//
//  Convert a (double) Modified Julian Day to GPS week and seconds of week.
//  The function returns a fractional gps week
// ---------------------------------------------------------------------------
function mjd_to_gps(dmjd)
{
    var mjd  = Math.floor(dmjd);
    var fmjd = dmjd - Math.floor(dmjd);
    var gps_week = (mjd - 44244)/7;
    var sec_of_week = ((mjd-44244)-gps_week*7+fmjd)*86400;
    return gps_week + sec_of_week / (7*86400);
}

// ---------------------------------------------------------------------------
//  function mjd_to_ymdhms
//
//  Convert a (double) Modified Julian Day to a string of type:
//  YYYY-MM-DD HH:MM:SS
// ---------------------------------------------------------------------------
function mjd_to_ymdhms(dmjd)
{
    var mjd  = Math.floor(dmjd);
    var fmjd = dmjd - Math.floor(dmjd);
    var month_day = [
        [0, 0],
        [31, 31],
        [59, 60],
        [90, 91],
        [120, 121],
        [151, 152],
        [181, 182],
        [212, 213],
        [243, 244],
        [273, 274],
        [304, 305],
        [334, 335],
        [365, 366]
    ];
    var days_fr_jan1_1901 = mjd - 15385;
    var num_four_yrs = parseInt(days_fr_jan1_1901/1461);
    var years_so_far = 1901 + parseInt(4*num_four_yrs);
    var days_left = days_fr_jan1_1901 - 1461*num_four_yrs;
    var delta_yrs = parseInt(days_left/365) - parseInt(days_left/1460);

    var year = years_so_far + delta_yrs;
    var yday = days_left - 365*delta_yrs + 1;
    var hour = parseInt(fmjd*24.0);
    var minute = parseInt(fmjd*1440.0 - hour*60.0);
    var second = parseInt(fmjd*86400.0 - hour*3600.0 - minute*60.0);
    var leap = 0;
    if (parseInt(year%4) == 0) leap = 1;
    var guess = parseInt(yday*0.032);
    var more = (( yday - month_day[guess+1][leap] ) > 0);
    var month = guess + more + 1;
    var mday = yday - month_day[guess+more][leap];

    return parseInt(year).toString() +
        "-" + parseInt(month).toString() +
        "-" + parseInt(mday).toString() +
        " " + parseInt(hour).toString() +
        ":" + parseInt(minute).toString() +
        ":" + parseInt(second).toString();
}

//  ---------------------------------------------------------------------------
//  function plot_raw()
//
//  This function will:
//  1. Parse the input json file (i.e. the time-series file and load data (locally)
//  2. Reset the x and y axis ranges and scales
//  3. Plot data as circles in three seperate SVGs; outliers are filled with
//+    yellow color.
//  4. Assign global min/max mjd
//  5. Enable/Disable Outliers and Model buttons depending on whether the model
//+    and events files are available (via globals g_md_filename and g_ev_filename)
//
// ----------------------------------------------------------------------------
function plot_raw()
{
    console.log("[DEBUG] Trying to parse and plot (raw) file \"" + g_ts_filename + "\".");
    
    //  Load data from file
    d3.json(g_ts_filename, function(error, data) {
        if ( error ) {
            window.alert("Failed to parse JSON file \""+ g_ts_filename + "\"");
            console.warn(error);
            throw error;
        }
        
        //  The actual data to plot is:
        pdata = data.data;

        //  Set x and y domains
        x.domain(d3.extent (pdata, function(d) { return d.epoch; }));
        yN.domain([0, 0]);
        yN.domain(d3.extent(pdata, function(d) { return d.north; }));
        yE.domain([0, 0]);
        yE.domain(d3.extent(pdata, function(d) { return d.east; }));
        yU.domain([0, 0]);
        yU.domain(d3.extent(pdata, function(d) { return d.up; }));
        
        //  Rescale y-axis
        svgN.selectAll(".y.axis").call(yAxisN);
        svgE.selectAll(".y.axis").call(yAxisE);
        svgU.selectAll(".y.axis").call(yAxisU);
        
        /*  Remove outliers if any
        svgN.selectAll(".residual_pt").remove();
        svgE.selectAll(".residual_pt").remove();
        svgU.selectAll(".residual_pt").remove();
        */
    
        //  Add the data points to all three svg (n, e, u).
        //  Points are added as circles; the normal color is red, though if a
        //+ data point is marked as outlier, it is filled with yellow color.
        //  All data points are marked as belonging to the class 'data_pt'. If
        //+ the data point is marked as outlier, it also belongs to the class
        //+ 'outlier_pt', else it also belongs to the class 'normal_pt'.
        svgN.selectAll(".dot")
            .data(pdata)
          .enter().append("circle")
            .attr("r", 3.0)
            .attr("cx", function(d) {return x(d.epoch);})
            .attr("cy", function(d) {return yN(d.north);})
            .attr("stroke", "blue")
            .attr("stroke-width", 0.3)
            .attr("fill", function(d) {
                if (d.flag_north.indexOf("o")>-1) {
                    return "yellow";
                }
                return "red";
            })
            .attr("class", function(d) {
                if (d.flag_north.indexOf("o")>-1) {
                    return "outlier_pt data_pt";
                } else {
                    return "normal_pt data_pt";
                }
            });
        svgE.selectAll(".dot")
            .data(pdata)
          .enter().append("circle")
            .attr("r", 3.0)
            .attr("cx", function(d) {return x(d.epoch);})
            .attr("cy", function(d) {return yN(d.east);})
            .attr("stroke", "blue")
            .attr("stroke-width", 0.3)
            .attr("fill", function(d) {
                if (d.flag_east.indexOf("o")>-1) {
                    return "yellow";
                }
                return "red";
            })
            .attr("class", function(d) {
                if (d.flag_east.indexOf("o")>-1) {
                    return "outlier_pt data_pt";
                } else {
                    return "normal_pt data_pt";
                }
            });
        svgU.selectAll(".dot")
            .data(pdata)
          .enter().append("circle")
            .attr("r", 3.0)
            .attr("cx", function(d) {return x(d.epoch);})
            .attr("cy", function(d) {return yN(d.up);})
            .attr("stroke", "blue")
            .attr("stroke-width", 0.3)
            .attr("fill", function(d) {
                if (d.flag_up.indexOf("o")>-1) {
                    return "yellow";
                }
                return "red";
            })
        .attr("class", function(d) {
            if (d.flag_up.indexOf("o")>-1) {
                return "outlier_pt data_pt";
            } else {
                return "normal_pt data_pt";
            }
        });
        
        //  Assign global min/max mjd
        g_min_mjd = x.domain()[0];
        g_max_mjd = x.domain()[1];
    });

    //  Enable/disable buttons
    document.getElementById("outlier-btn").disabled = false;
    if ( g_md_filename.trim() ) {
        document.getElementById("model-btn").disabled   = false;
    } else {
        document.getElementById("model-btn").disabled   = true;
    }
    if ( g_ev_filename.trim() ) {
        document.getElementById("events-btn").disabled  = false;
    } else {
        document.getElementById("events-btn").disabled  = true;
    }
    
    console.log("[DEBUG] File \"" + g_ts_filename + "\" parsed and plotted.");
}

// ----------------------------------------------------------------------------
//  function remove_datapt
//
//  Remove all data points belonging to the class 'data_pt' from all three SVGs
//+ (svgN, svgE, svgU)
//
// ----------------------------------------------------------------------------
function remove_datapt()
{
    svgN.selectAll(".data_pt").remove();
    svgE.selectAll(".data_pt").remove();
    svgU.selectAll(".data_pt").remove();
}

// ----------------------------------------------------------------------------
//  function plot_residuals
//
//
//  Read in a json time-series file and plot residuals (North, East and Up).
//  The function will remove any points/elements of the plots, belonging to
//+ the class 'data_pt'. All points ploted (on the svg's) will be of class
//+ 'residual_pt'
//  The plots will be re-scaled and outliers will be filtered out (they will actually
//+ be plotted as circles with radius zero).
//  The y-axis will be set to MilliMeters.
//  This function will also set the 'Filter Outliers' button to disabled state.
//  TODO Rename the y-axis
//
// ----------------------------------------------------------------------------
function plot_residuals()
{
    console.log("[DEBUG] Plotting residuals using file \"" + g_ts_filename + "\".");
    
    d3.json(g_ts_filename, function(error, data) {
        if ( error ) throw error;
        pdata = data.data;
        
        //  Set x and y domains
        x.domain(d3.extent (pdata, function(d) { return d.epoch; }));
        yN.domain([0, 0]);
        yN.domain(d3.extent(pdata, function(d) {
            return (d.flag_north.indexOf("o")>-1) ? 0 : d.res_north*1000;
        }));
        yE.domain([0, 0]);
        yE.domain(d3.extent(pdata, function(d) {
            return (d.flag_east.indexOf("o")>-1) ? 0 : d.res_east*1000;
        }));
        yU.domain([0, 0]);
        yU.domain(d3.extent(pdata, function(d) {
            return (d.flag_up.indexOf("o")>-1) ? 0 : d.res_up*1000;
        }));
        
        //  Rescale y-axis
        svgN.selectAll(".y.axis").call(yAxisN);
        svgE.selectAll(".y.axis").call(yAxisE);
        svgU.selectAll(".y.axis").call(yAxisU);

        //  Remove any data points belonging to the class 'data_pt'
        svgN.selectAll(".data_pt").remove();
        svgE.selectAll(".data_pt").remove();
        svgU.selectAll(".data_pt").remove();
        
        //  Plot the residuals as points (circles)
        svgN.selectAll(".dot")
            .data(pdata)
          .enter().append("circle")
            .attr("r", function(d) {
                return (d.flag_north.indexOf("o")>-1) ? 0 : 3.0;
            })
            .attr("cx", function(d) {return x(d.epoch);})
            .attr("cy", function(d) {return yN(d.res_north * 1000);})
            .attr("stroke", "blue")
            .attr("stroke-width", 0.3)
            .attr("fill", "red")
            .attr("class", function(d) {
                return (d.flag_north.indexOf("o")>-1)
                    ? "outlier_pt residual_pt"
                    : "normal_pt residual_pt";
            });
        svgE.selectAll(".dot")
            .data(pdata)
          .enter().append("circle")
            .attr("r", function(d) {
                return (d.flag_east.indexOf("o")>-1) ? 0 : 3.0;
            })
            .attr("cx", function(d) {return x(d.epoch);})
            .attr("cy", function(d) {return yE(d.res_east * 1000);})
            .attr("stroke", "blue")
            .attr("stroke-width", 0.3)
            .attr("fill", "red")
            .attr("class", function(d) {
                return (d.flag_east.indexOf("o")>-1)
                    ? "outlier_pt residual_pt"
                    : "normal_pt residual_pt";
            });
        svgU.selectAll(".dot")
            .data(pdata)
          .enter().append("circle")
            .attr("r", function(d) {
                return (d.flag_up.indexOf("o")>-1) ? 0 : 3.0;
            })
            .attr("cx", function(d) {return x(d.epoch);})
            .attr("cy", function(d) {return yU(d.res_up * 1000);})
            .attr("stroke", "blue")
            .attr("stroke-width", 0.3)
            .attr("fill", "red")
            .attr("class", function(d) {
                return (d.flag_up.indexOf("o")>-1)
                    ? "outlier_pt residual_pt"
                    : "normal_pt residual_pt";
            });
        
    });
    
    //  Check/uncheck buttons button
    document.getElementById("outlier-btn").checked = false;
    document.getElementById("model-btn").checked   = false;
    
    //  Disable outlier button
    document.getElementById("outlier-btn").disabled = true;
    document.getElementById("model-btn").disabled = true;
    
    console.log("[DEBUG] Residuals plotted");
}

// ---------------------------------------------------------------------------
//  function plot_raw_no_outliers
//
//  This function will:
//  1. Remove all points (on all three svg's) tha belong to the 'outlier_pt'
//+    class
//  2. Rescale the y-axis (to not include the above points)
//  3. Re-scale and re-plot all points belonging to the 'normal_pt' class
//
// ---------------------------------------------------------------------------
function plot_raw_no_outliers()
{
    console.log("[DEBUG] plotting time-series with no outliers, file is \""+g_ts_filename+"\".");

    if ( document.getElementById("outlier-btn").checked === false ) {
        svgN.selectAll(".data_pt").remove();
        svgE.selectAll(".data_pt").remove();
        svgU.selectAll(".data_pt").remove();
        plot_raw();
        return;
    }

    // Get the data again
    d3.json(g_ts_filename, function(error, data) {
        if ( error ) throw error;
        pdata = data.data;

        // set x and y domains
        x.domain(d3.extent (pdata, function(d) { return d.epoch; }));
        yN.domain([0, 0]);
        yE.domain([0, 0]);
        yU.domain([0, 0]);
        yN.domain(d3.extent(pdata, function(d) {
            return (d.flag_north.indexOf("o")>-1) ? 0 : d.north;
        }));
        yE.domain(d3.extent(pdata, function(d) { 
            return (d.flag_east.indexOf("o")>-1) ? 0 : d.east;
        }));
        yU.domain(d3.extent(pdata, function(d) { 
            return (d.flag_up.indexOf("o")>-1) ? 0 : d.up;
        }));
        
        // rescale y-axis
        svgN.selectAll(".y.axis").call(yAxisN);
        svgE.selectAll(".y.axis").call(yAxisE);
        svgU.selectAll(".y.axis").call(yAxisU);
    
        // remove points that belong to the 'outlier_pt' class
        svgN.selectAll(".outlier_pt").remove();
        svgE.selectAll(".outlier_pt").remove();
        svgU.selectAll(".outlier_pt").remove();
        
        // re-scale points that belong to the 'normal_pt' class
        svgN.selectAll(".data_pt")
            .attr("cx", function(d) {return x(d.epoch);})
            .attr("cy", function(d) {return yN(d.north);});
        svgE.selectAll(".data_pt")
            .attr("cx", function(d) {return x(d.epoch);})
            .attr("cy", function(d) {return yE(d.east);});
        svgU.selectAll(".data_pt")
            .attr("cx", function(d) {return x(d.epoch);})
            .attr("cy", function(d) {return yU(d.up);});
            
        // if we also have the model lines, we need to rescale them!
        if ( svgN.select(".model_line").empty() ) {
            console.log("empty selection");
        } else {
            console.log("rescaling model lines");
            // store model points here
            var xdata, ydata, zdata;
            // Get the data
            d3.json(g_md_filename, function(error, data) {
                if ( error ) throw error;
                
                // compute the model points/canvas
                xdata = make_model(data.model_x, g_min_mjd, g_max_mjd, 1);
                ydata = make_model(data.model_y, g_min_mjd, g_max_mjd, 1);
                zdata = make_model(data.model_z, g_min_mjd, g_max_mjd, 1);

                // re-scale the model-line (if any)
                svgN.select(".model_line").attr("d", nline(xdata));
                svgE.select(".model_line").attr("d", eline(ydata));
                svgU.select(".model_line").attr("d", uline(zdata));
            });
        }
    });
}

// ---------------------------------------------------------------------------
//  function make_model
//
//  Given a "model" object (as read off from a corrsponding json file),
//+ compute the values of the model for values in range [from, to] with a 
//+ step = step.
//  The values are returned as an array of objects of type:
//+ [...{"t":...,"val":...}...] 
//+ so that they can be ploted via the [neu]line functions.
//
//// -------------------------------------------------------------------------
function make_model(mod, from, to, step)
{
    console.log("[DEBUG] Making model line ...");

    var data = [];
    var value = 0;
    var ref_t = mod.reference_epoch;

    for (var t = from; t <= to; t += step) {
        var dt = t - ref_t;
        value  = mod.const_term;
        value += mod.velocity*dt/365.25;
        for (var comp in mod.harmonics) {
            if (t>=mod.harmonics[comp].from && t<mod.harmonics[comp].to)
            {
                var angular_frequency = 2*Math.PI/mod.harmonics[comp].period;
                value += mod.harmonics[comp].in_phase*Math.cos(angular_frequency*dt);
                value += mod.harmonics[comp].out_of_phase*Math.sin(angular_frequency*dt);
            }
        }
        for (var comp in mod.jumps) {
            if (t>=mod.jumps[comp].at) {
                value += mod.jumps[comp].value;
            }
        }
        data.push({"t":t, "val":value});
    }
    return data;
}

function remove_models()
{
    svgN.selectAll(".model_line").remove();
    svgE.selectAll(".model_line").remove();
    svgU.selectAll(".model_line").remove();
}
/*
** Read coordinate time-series model(s) off from a corresponding json file,
** and plot the models (one per component). The plotted lines will belong to a
** class(es) named "line model_line"
*/
function plot_models()
{
    console.log("[DEBUG] Plotting model line using file \""+g_md_filename+"\".");

    if ( document.getElementById("model-btn").checked === false ) {
        remove_models();
        return;
    }

    // store model points here
    var xdata, ydata, zdata;
    
    // Get the data
    d3.json(g_md_filename, function(error, data) {
        if ( error ) throw error;
        
        // compute the model points/canvas
        xdata = make_model(data.model_x, g_min_mjd, g_max_mjd, 1);
        ydata = make_model(data.model_y, g_min_mjd, g_max_mjd, 1);
        zdata = make_model(data.model_z, g_min_mjd, g_max_mjd, 1);
        
        // make model a list of strings
        xmodel = data.model_x;
        ymodel = data.model_y;
        zmodel = data.model_z;
        xstr = model2text(xmodel);
        ystr = model2text(ymodel);
        zstr = model2text(zmodel);

        // Define the div for the tooltip, see http://bl.ocks.org/d3noob/a22c42db65eb00d4e369
        var tooltip_div = d3.select("body").append("div")   
                            .attr("class", "model-tooltip")               
                            .style("opacity", 0);

        // append the model as line
        svgN.append("path")
            .datum(xdata)
            .attr("class", "line model_line")
            .attr("d", nline)
            .on("mouseover", function(d) {
                tooltip_div.transition()
                           .duration(200)
                           .style("opacity", .9);
                tooltip_div.html(xstr)
                           .style("left", (d3.event.pageX) + "px")
                           .style("top", (d3.event.pageY - 28) + "px");
            })
            .on("mouseout", function(d) {
                tooltip_div.transition()
                           .duration(200)
                           .style("opacity", 0);
            });
        svgE.append("path")
            .datum(ydata)
            .attr("class", "line model_line")
            .attr("d", eline)
            .on("mouseover", function(d) {
                tooltip_div.transition()
                           .duration(200)
                           .style("opacity", .9);
                tooltip_div.html(ystr)
                           .style("left", (d3.event.pageX) + "px")
                           .style("top", (d3.event.pageY - 28) + "px");
            })
            .on("mouseout", function(d) {
                tooltip_div.transition()
                           .duration(200)
                           .style("opacity", 0);
            });
        svgU.append("path")
            .datum(zdata)
            .attr("class", "line model_line")
            .attr("d", uline)
            .on("mouseover", function(d) {
                tooltip_div.transition()
                           .duration(200)
                           .style("opacity", .9);
                tooltip_div.html(zstr)
                           .style("left", (d3.event.pageX) + "px")
                           .style("top", (d3.event.pageY - 28) + "px");
            })
            .on("mouseout", function(d) {
                tooltip_div.transition()
                           .duration(200)
                           .style("opacity", 0);
            });
    });
}

function model2text(model)
{
    /* As string with newlines
       ------------------------
    var str = 'Reference Epoch: ';
    str += mjd_to_ymdhms(model.reference_epoch) + '\n';
    str += 'Velocity       : ';
    str += model.velocity + ' m/yr\n';
    for (var i in model.harmonics) {
        str += 'Period         : ' + model.harmonics[i].period;
        str += '\n  In Phase     : ' + model.harmonics[i].in_phase + ' (m)';
        str += '\n  Out Of Phase : ' + model.harmonics[i].out_of_phase + ' (m)';
    }
    for (var i in model.jumps) {
        str += '\nJump           : ' + model.jumps[i].value + ' (m) at ' + mjd_to_ymdhms(model.jumps[i].at);
    }
    return str;
    */
    /* As list of strings (lines)
       ----------------------------
    var line_list = [];
    line_list.push('Reference Epoch: ' + mjd_to_ymdhms(model.reference_epoch) + '\n');
    line_list.push('Velocity       : ' + model.velocity + ' m/yr\n');
    for (var i in model.harmonics) {
        line_list.push('Period         : ' + model.harmonics[i].period);
        line_list.push('  In Phase     : ' + model.harmonics[i].in_phase + ' (m)');
        line_list.push('  Out Of Phase : ' + model.harmonics[i].out_of_phase + ' (m)');
    }
    for (var i in model.jumps) {
        line_list.push('Jump           : ' + model.jumps[i].value + ' (m) at ' + mjd_to_ymdhms(model.jumps[i].at));
    }
    return line_list;
    */
    var str = "";
    str += '<p>Reference Epoch: ';
    str += mjd_to_ymdhms(model.reference_epoch) + '</p>';
    str += '<p>Velocity       : ';
    str += model.velocity + ' m/yr</p>';
    for (var i in model.harmonics) {
        str += '<p>Period         : ' + model.harmonics[i].period;
        str += '<br>  In Phase     : ' + model.harmonics[i].in_phase + ' (m)';
        str += '<br>  Out Of Phase : ' + model.harmonics[i].out_of_phase + ' (m)</p>';
    }
    for (var i in model.jumps) {
        str += '<p>Jump           : ' + model.jumps[i].value + ' (m) at ' + mjd_to_ymdhms(model.jumps[i].at)+"</p>";
    }
    return str;
}
/*
function show_model_info()
{
    // Get the data
    d3.json(g_md_filename, function(error, data) {
        if ( error ) throw error;

        xmodel = data.model_x;
        ymodel = data.model_y;
        zmodel = data.model_z;

        xstr = model2text(xmodel);
        var start_x, start_y;
        for (var i in xstr) {
            svgN.append("text")
                .attr("class", "model_text")
                .attr("x", x(start_x))
                .attr("y", yN(start_y))
                .text(xstr[i]);
            start_y -= 0.001;
        }
        //console.log(xstr);
    });
}
*/

function remove_events()
{
    svgN.selectAll(".event_line").remove();
    svgE.selectAll(".event_line").remove();
    svgU.selectAll(".event_line").remove();
}
function plot_events()
{
    console.log("[DEBUG] Plotting events using file \""+g_ev_filename+"\".");
    if ( document.getElementById("events-btn").checked === false ) {
        remove_events();
        return;
    }

    // Get the data
    d3.json(g_ev_filename, function(error, data) {
        if ( error ) throw error;
        
        // split events
        jump_list   = data.jumps;
        velchg_list = data.velocity_changes;
        erthqk_list = data.earthquakes;

        for (var evnt in jump_list) {
            var t  = jump_list[evnt].at;
            svgN.append("line")
                .attr("x1", x(t))
                .attr("y1", yN(yAxisN.scale().domain()[0]))
                .attr("x2", x(t))
                .attr("y2", yN(yAxisN.scale().domain()[1]))
                .attr("class", "event_line jump_line");
            svgE.append("line")
                .attr("x1", x(t))
                .attr("y1", yE(yAxisE.scale().domain()[0]))
                .attr("x2", x(t))
                .attr("y2", yE(yAxisE.scale().domain()[1]))
                .attr("class", "event_line jump_line");
            svgU.append("line")
                .attr("x1", x(t))
                .attr("y1", yU(yAxisU.scale().domain()[0]))
                .attr("x2", x(t))
                .attr("y2", yU(yAxisU.scale().domain()[1]))
                .attr("class", "event_line jump_line");
        }
        
        for (var evnt in velchg_list) {
            var t  = velchg_list[evnt].at;
            svgN.append("line")
                .attr("x1", x(t))
                .attr("y1", yN(yAxisN.scale().domain()[0]))
                .attr("x2", x(t))
                .attr("y2", yN(yAxisN.scale().domain()[1]))
                .attr("class", "event_line velchg_line");
            svgE.append("line")
                .attr("x1", x(t))
                .attr("y1", yE(yAxisE.scale().domain()[0]))
                .attr("x2", x(t))
                .attr("y2", yE(yAxisE.scale().domain()[1]))
                .attr("class", "event_line velchg_line");
            svgU.append("line")
                .attr("x1", x(t))
                .attr("y1", yU(yAxisU.scale().domain()[0]))
                .attr("x2", x(t))
                .attr("y2", yU(yAxisU.scale().domain()[1]))
                .attr("class", "event_line velchg_line");
        }
        for (var evnt in erthqk_list) {
            var t  = erthqk_list[evnt].at;
            svgN.append("line")
                .attr("x1", x(t))
                .attr("y1", yN(yAxisN.scale().domain()[0]))
                .attr("x2", x(t))
                .attr("y2", yN(yAxisN.scale().domain()[1]))
                .attr("class", "event_line erthqk_line");
            svgE.append("line")
                .attr("x1", x(t))
                .attr("y1", yE(yAxisE.scale().domain()[0]))
                .attr("x2", x(t))
                .attr("y2", yE(yAxisE.scale().domain()[1]))
                .attr("class", "event_line erthqk_line");
            svgU.append("line")
                .attr("x1", x(t))
                .attr("y1", yU(yAxisU.scale().domain()[0]))
                .attr("x2", x(t))
                .attr("y2", yU(yAxisU.scale().domain()[1]))
                .attr("class", "event_line erthqk_line");
        }
    });
}

function plot_residual_histogram()
{
    console.log("[DEBUG] Plotting residual histogram using file \""+g_ts_filename+"\".");
    remove_events();
    remove_models();
    remove_datapt();


    d3.json(g_ts_filename, function(error, data) {
        if ( error ) throw error;
        pdata = data.data;
     
        var xN = d3.scaleLinear().range([0, width]);
        var xAxisN  = d3.axisBottom(xN);
        var xE = d3.scaleLinear().range([0, width]);
        var xAxisE  = d3.axisBottom(xE);
        var xU = d3.scaleLinear().range([0, width]);
        var xAxisU  = d3.axisBottom(xU);
        var x_axis_limits = {"xnmin":0, "xnmax":0, "xemin":0, "xemax":0, "xumin":0, "xumax":0};
        for (var o in data.data) {
            if (data.data[o].flag_north.indexOf("o") === -1) {
                if (x_axis_limits.xnmin > data.data[o].res_north)
                    x_axis_limits.xnmin = data.data[o].res_north;
                if (x_axis_limits.xnmax < data.data[o].res_north)
                    x_axis_limits.xnmax = data.data[o].res_north;
            }
            else {
                console.log("excluding outlier:", data.data[o].res_north);
            }
            if (data.data[o].flag_east.indexOf("o") === -1) {
                if (x_axis_limits.xemin > data.data[o].res_east)
                    x_axis_limits.xemin = data.data[o].res_east;
                if (x_axis_limits.xemax < data.data[o].res_east)
                    x_axis_limits.xemax = data.data[o].res_east;
            }
            if (data.data[o].flag_up.indexOf("o") === -1) {
                if (x_axis_limits.xumin > data.data[o].res_up)
                    x_axis_limits.xumin = data.data[o].res_up;
                if (x_axis_limits.xumax < data.data[o].res_up)
                    x_axis_limits.xumax = data.data[o].res_up;
            }
        }


        xN.domain([x_axis_limits.xnmin, x_axis_limits.xnmax]);
        xE.domain([x_axis_limits.xemin, x_axis_limits.xemax]);
        xU.domain([x_axis_limits.xumin, x_axis_limits.xumax]);
        
        var Nhistogram = d3.histogram()
                        .domain(xN.domain())
                        .thresholds(xN.ticks(20));
        var Ehistogram = d3.histogram()
                        .domain(xN.domain())
                        .thresholds(xN.ticks(20));
        var Uhistogram = d3.histogram()
                        .domain(xN.domain())
                        .thresholds(xN.ticks(20));
        
        ResN = data.data.map(function(a) {;
            return (a.flag_north.indexOf("o")>-1) ? a.res_north : null;
        });
        ResE = data.data.map(function(a) {
            return (a.flag_east.indexOf("o")>-1) ? a.res_east : null;
        });
        ResU = data.data.map(function(a) {
            return (a.flag_up.indexOf("o")>-1) ? a.res_up : null;
        });

        var BinN = Nhistogram(ResN);
        var BinE = Ehistogram(ResE);
        var BinU = Uhistogram(ResU);
        
        yN.domain([0, d3.max(BinN, function(d) { return d.length; })]);
        yE.domain([0, d3.max(BinE, function(d) { return d.length; })]);
        yU.domain([0, d3.max(BinU, function(d) { return d.length; })]);
        svgN.selectAll(".y.axis").call(yAxisN);
        svgE.selectAll(".y.axis").call(yAxisE);
        svgU.selectAll(".y.axis").call(yAxisU);
        svgN.selectAll(".x.axis").call(xAxisN);
        svgE.selectAll(".x.axis").call(xAxisE);
        svgU.selectAll(".x.axis").call(xAxisU);

        svgN.selectAll("rect")
                .data(BinN)
            .enter().append("rect")
                .attr("class", "bar")
                .attr("x", 1)
                .attr("transform", function(d) {
                    return "translate(" + xN(d.x0) + "," + yN(d.length) + ")"; })
                .attr("width", function(d) { return xN(d.x1) - xN(d.x0) -1 ; })
                .attr("height", function(d) { return height - yN(d.length); });
        svgE.selectAll("rect")
                .data(BinE)
            .enter().append("rect")
                .attr("class", "bar")
                .attr("x", 1)
                .attr("transform", function(d) {
                    return "translate(" + xE(d.x0) + "," + yE(d.length) + ")"; })
                .attr("width", function(d) { return xE(d.x1) - xE(d.x0) -1 ; })
                .attr("height", function(d) { return height - yE(d.length); });
        svgU.selectAll("rect")
                .data(BinU)
            .enter().append("rect")
                .attr("class", "bar")
                .attr("x", 1)
                .attr("transform", function(d) {
                    return "translate(" + xU(d.x0) + "," + yU(d.length) + ")"; })
                .attr("width", function(d) { return xU(d.x1) - xU(d.x0) -1 ; })
                .attr("height", function(d) { return height - yU(d.length); });
    });
}

/*
** Attach the x axis to a datetime scale in the domain [g_min_mjd,g_max_mjd]
** Note that nothing gets reploted!
*/
function to_ymdhms()
{
    x = d3.scaleTime().range([0, width]);
    min_ymd = mjd_to_ymdhms(g_min_mjd);
    max_ymd = mjd_to_ymdhms(g_max_mjd);
    // console.log("Min date:", min_ymd);
    // console.log("Max date:", max_ymd);
    x.domain([dateFormat.parse(min_ymd), dateFormat.parse(max_ymd)]);
    xAxis  = d3.axisBottom(x);
    
    // rescale x-axis
    svgN.selectAll(".x.axis").call(xAxis);
    svgE.selectAll(".x.axis").call(xAxis);
    svgU.selectAll(".x.axis").call(xAxis);
}
/*
** Attach the x axis to a datetime scale in the domain [g_min_mjd,g_max_mjd]
** Note that nothing gets reploted!
*/
function to_mjd()
{
    x = d3.scaleLinear().range([0, width]); // v4.x
    x.domain([g_min_mjd, g_max_mjd]);
    xAxis  = d3.axisBottom(x);
    
    // rescale x-axis
    svgN.selectAll(".x.axis").call(xAxis);
    svgE.selectAll(".x.axis").call(xAxis);
    svgU.selectAll(".x.axis").call(xAxis);
}
function to_gpsw()
{
    x = d3.scaleLinear().range([0, width]); // v4.x
    x.domain([mjd_to_gps(g_min_mjd), mjd_to_gps(g_max_mjd)]);
    xAxis  = d3.axisBottom(x);
    
    // rescale x-axis
    svgN.selectAll(".x.axis").call(xAxis);
    svgE.selectAll(".x.axis").call(xAxis);
    svgU.selectAll(".x.axis").call(xAxis);
}

// ---------------------------------------------------------------------------
//  function files_changed
//
//  This function is triggered by the SUBMIT button.
//  1. First collect the input files filenames.
//  2. Remove any plots/lines that already exist on page.
//  3. Plot the raw time-series
//
// ---------------------------------------------------------------------------
function files_changed()
{
    console.log("[DEBUG] In function files_changed()");
    g_ts_filename = document.getElementById("ts-file-selector").value;
    g_md_filename = document.getElementById("md-file-selector").value;
    g_ev_filename = document.getElementById("ev-file-selector").value;

    remove_events();
    remove_models();
    remove_datapt();

    try {
        plot_raw();
    } catch(err) {
        console.log("OOPS!!");
        console.log("Got g_ts_filename = ", g_ts_filename);
        console.log("Got g_md_filename = ", g_md_filename);
        if ( !g_md_filename.trim() ) {
            console.log("empty model filename");
        } else {
            console.log("model file is not empty, but \""+g_md_filename.trim()+"\".");
        }
        console.log("Got g_ev_filename = ", g_ev_filename);
        if ( !g_ev_filename.trim() ) {console.log("empty event filename");}
        // window.alert("ERROR. Something went wrong in plotting the (raw) time-series.");
    }
    console.log("[DEBUG] Out function files_changed()");
}
