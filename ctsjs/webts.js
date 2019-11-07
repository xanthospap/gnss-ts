var DEBUG_MODE = true;

// ---------------------------------------------------------------------------
//
//  function remove_axis_labels
//
//  Remove axis lables on all three SVG's. Actually, the function will remove
//+ everything belonging o the class 'axis-label', attached to any of the three
//+ SVG.
//
// ---------------------------------------------------------------------------
function remove_axis_labels()
{
  svgN.selectAll(".axis-label").remove();
  svgE.selectAll(".axis-label").remove();
  svgU.selectAll(".axis-label").remove();
}

// ---------------------------------------------------------------------------
//  function add_axis_label
//
//  Add x- and y- axis labels to the three SVG's.
//  See http://bl.ocks.org/phoebebright/3061203
//
//  --------------------------------------------------------------------------
function add_axis_label(xlabel, yn_label, ye_label, yu_label)
{
  if (DEBUG_MODE) console.log("[DEBUG] fun::add_axis_label() Igniting ...");
  // space around the chart, not including labels
  padding = -1.5 * margin.left;

  //  First remove any already-in-place axis labels
  remove_axis_labels();
  
  svgN.append("text")
  .attr("text-anchor", "middle")
  .attr("transform", "translate("+ (width/2) +","+(height-(padding/2))+")")
  .attr("class", "axis-label")
  .text(xlabel);
  svgN.append("text")
  .attr("text-anchor", "middle")
  .attr("transform", "translate("+ (padding/2) +","+(height/2)+")rotate(-90)")
  .attr("class", "axis-label")
  .text(yn_label);
  
  svgE.append("text")
  .attr("text-anchor", "middle")
  .attr("transform", "translate("+ (width/2) +","+(height-(padding/2))+")")
  .attr("class", "axis-label")
  .text(xlabel);
  svgE.append("text")
  .attr("text-anchor", "middle")
  .attr("transform", "translate("+ (padding/2) +","+(height/2)+")rotate(-90)")
  .attr("class", "axis-label")
  .text(ye_label);
  
  svgU.append("text")
  .attr("text-anchor", "middle")
  .attr("transform", "translate("+ (width/2) +","+(height-(padding/2))+")")
  .attr("class", "axis-label")
  .text(xlabel);
  svgU.append("text")
  .attr("text-anchor", "middle")
  .attr("transform", "translate("+ (padding/2) +","+(height/2)+")rotate(-90)")
  .attr("class", "axis-label")
  .text(yu_label);
  
  //  Set global axis names
  g_xlabel   = xlabel;
  g_yn_label = yn_label;
  g_ye_label = ye_label;
  g_yu_label = yu_label;
  
  if (DEBUG_MODE) console.log("[DEBUG] fun::add_axis_label() All done");
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
  if (DEBUG_MODE) console.log("[DEBUG] fun::plot_raw() Igniting ...");
  
  //  If we are showing the residual plot, we need to remove it first
  // remove_residuals();
  
  //  Load data from file
  if (DEBUG_MODE) console.log("\t[DEBUG] fun::plot_raw() Loading/parsing file: ", g_ts_filename);
  d3.json(g_ts_filename, function(error, data) {
    if ( error ) throw error;
          
    //  The actual data to plot is:
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async -- ", data[0]);
    pdata = data;
    
    //  Set x and y domains
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async --  Setting x, N, E, U domains");
    x.domain(d3.extent (pdata, function(d) { return parseFloat(d.epoch); }));
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async --  New x domain", x.domain());
    yN.domain([0, 0]);
    yN.domain(d3.extent(pdata, function(d) { return parseFloat(d.north); }));
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async --  New N domain", yN.domain());
    yE.domain([0, 0]);
    yE.domain(d3.extent(pdata, function(d) { return parseFloat(d.east); }));
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async --  New E domain", yE.domain());
    yU.domain([0, 0]);
    yU.domain(d3.extent(pdata, function(d) { return parseFloat(d.up); }));
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async --  New U domain", yU.domain());
          
    //  Rescale y-axis
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async --  Rescaling axis");
    svgN.select(".x.axis").call(xAxis);
    svgN.selectAll(".y.axis").call(yAxisN);
    svgE.selectAll(".y.axis").call(yAxisE);
    svgU.selectAll(".y.axis").call(yAxisU);
    
    //  Add axis labels
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async --  Adding axis labels");
    add_axis_label( g_xlabel, "North (m)", "East (m)", "Up (m)" );
    
    //  Add the data points to all three svg (n, e, u).
    //  Points are added as circles; the normal color is red, though if a
    //+ data point is marked as outlier, it is filled with yellow color.
    //  All data points are marked as belonging to the class 'data_pt'. If
    //+ the data point is marked as outlier, it also belongs to the class
    //+ 'outlier_pt', else it also belongs to the class 'normal_pt'.
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async -- Plotting N component");
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
          
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async -- Plotting E component");
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
          
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async -- Plotting U component");
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
    if (DEBUG_MODE) console.log("\t[DEBUG] fun:plot_raw() -- async -- Assigning global min/max mjd");
    g_min_mjd = x.domain()[0];
    g_max_mjd = x.domain()[1];
  });
  
  //  Enable/disable buttons
  if (DEBUG_MODE) console.log("\t[DEBUG] fun::plot_raw() Enable/Disable buttons");
  /*
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
  }*/
  if (DEBUG_MODE) console.log("[DEBUG] fun::plot_raw() All done");
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
  if (DEBUG_MODE) console.log("[DEBUG] fun::files_changed() Ignited...");
  g_ts_filename = document.getElementById("time-series-file-selector").value;
  //g_md_filename = document.getElementById("md-file-selector").value;
  //g_ev_filename = document.getElementById("ev-file-selector").value;
  g_ts_filename = "output.json"
  
  /*
  remove_events();
  remove_models();
  remove_datapt();
  */
  
  try {
    plot_raw();
  } catch(err) {
    console.log("\t[DEBUG] fun::files_changed() Error caught!");
  }
  console.log("[DEBUG] fun::files_changed() All done");
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
          pdata = data;
    
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
