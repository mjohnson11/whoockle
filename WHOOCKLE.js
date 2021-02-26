var whoo_holder;

var whoo_DEBUG = true;

var whoo_colors = {
  'Positive': 'hsl(350, 100%, 33%)',
  'Negative': 'hsl(0, 0%, 33%)',
  'Inconclusive': 'hsl(34, 100%, 33%)',
  'Rerun-A': 'hsl(34, 100%, 33%)',
  'Rerun-B': 'hsl(34, 100%, 33%)'
};

var whoo_gene_colors = {
  'N1': '#618BE0', // reporter is FAM which is blue/green
  'RDRP': '#00E66C',   // reporter is VIC which is green
  'RNASEP': '#FF2416'  // reporter is CY5 which is red
};

var whoo_export_map = {
  'Inconclusive': 'INCONCLUSIVE RESULT for SARS-CoV-2',
  'Negative': 'NEGATIVE for SARS-CoV-2',
  'Positive': 'POSITIVE for SARS-CoV-2',
  'Rerun-A': 'Retest A',
  'Rerun-B': 'Retest B'
}

var control_names = ['PTC', 'NTC', 'NTC-W', 'NTC-P']; // possible sample names for controls

var whoo_table_row_height = 20;

var whoo_possible_results = ['Positive', 'Inconclusive', 'Rerun-A', 'Rerun-B', 'Negative'];

var whoo_rows = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P'];

class whooData {
  /**
   * This class holds all the data and variables associated with a single open file
   */
  constructor() {
  }
  
  //#region FILE/DATA HANDLING FUNCTIONS
  
  load_from_save() {
    /**
     * A function for taking already processed data and re-loading it
     */
  }

  save() {
    /**
     * A function for saving data in a format that can be reopened
     */
  }

  format_for_export() {
    /**
     * Formats a new dataframe for export
     * Columns: Well, Well Position, N1 CT, RDRP CT, RNASEP CT, RawCall, OverrideCall, Override, FinalCall
     */
    this.export_data = [];
    let notCleared = ['Rerun-A', 'Rerun-B'];
    //let simple_columns = ['RawCall', 'OverrideCall', 'Override'];
    for (let d of this.main_data) {
      let tmp_row = {'Position': d['Well Position'], 'InterpretiveResult': whoo_export_map[d['FinalCall']]};
      for (let c of ['N1 CT', 'RDRP CT', 'RNASEP CT']) {
        tmp_row[c.split(' ')[0]] = (d[c]) ? d[c].toFixed(1) : 'Undetermined'; // yeilds Undetermined if CT is NaN, the number formatted to one dec. point if not
      }
      // Only Positive, Negative, and Inconclusive calls, excluding control wells, are cleared to report
      tmp_row['ClearToReport'] = ( (notCleared.indexOf(d['FinalCall']) == -1) && (control_names.indexOf(d['Sample Name']) == -1) );
      //for (let c of simple_columns) { // adding a few extra columns, currently removed so parser can work easier
      //  tmp_row[c] = d[c];
      //}
      this.export_data.push(tmp_row);
    }
  }

  export() {
    /**
     * Exports final calls and CTs of each well
     */
    this.format_for_export();
    let output_name = 'Processed_' + this.input_file.name.replace('.txt', '.csv');
    let plate_ID = this.input_file.name.split(' ')[1].split('.')[0];
    let a = document.createElement('a');
    let export_header = 'QpcrPlateId: ' + plate_ID + '\n';
    export_header += 'InstrumentSerialNumber: ' + this.raw_data.header['Instrument Serial Number'] + '\n\n';
    let output_file = new Blob([export_header + d3.csvFormat(this.export_data)], {type: 'text/plain'});
    
    a.href= URL.createObjectURL(output_file);
    a.download = output_name;
    a.click();

    URL.revokeObjectURL(a.href);
  }

  load_file() {
    /**
     * Loads a QuantStudio 7 Flex txt output file
     * (see https://assets.thermofisher.com/TFS-Assets/LSG/manuals/4489825.pdf pgs. 59-66),
     * this function builds the raw_data object with properties for each part of the file:
     *  - One for the header information (stored as plain text)
     *  - One for each of the additional sections inside the file structure (parsed by d3.tsvParse)
     *  - For our use, these are "Sample Setup", "Amplification Data", and "Results"
     *  - we ignore anything appearing after a blank line in any of these sections
     *    (eg at the end of Results there is "Analysis Type" and "Reference Sample" info)
     */
    if (whoo_DEBUG) console.log('in load_file()');
    let fr = new FileReader(); 
    let self = this; // to avoid "this" pointing to the wrong thing inside the onload
    fr.onload=function() { // this function runs once the file is loaded - all the action is in here
      let file_string = fr.result;
      // a legacy replacement for old files where N1 was N GENE
      file_string = file_string.replaceAll('N GENE', 'N1')
      let split_by_blank_lines = file_string.split(/\n\s*\n/); // split file by blank lines
      self.raw_data = {'header': {}}; // Parsing header info...
      for (let line of split_by_blank_lines[0].split('\n')) {
        self.raw_data.header[line.split(' = ')[0].slice(2,)] = line.split(' = ')[1];
      }
      for (let block of split_by_blank_lines) {
        if (block[0]=='[') {
          let block_name = block.slice(1,block.indexOf(']'));
          console.log(block_name);
          self.raw_data[block_name] = d3.tsvParse(block.slice(block.indexOf('\n')+1));
        }
      }
      if (whoo_DEBUG) console.log('file loaded');
      self.make_main_data();
    }
    this.input_file = document.getElementById('inputfile').files[0];
    fr.readAsText(document.getElementById('inputfile').files[0]);
  }

  row_by_well(well) {
    /**
     * Returns the row of the dataframe associated with a specific well
     * "well" is the well number, though it is stored as a string
     */
    return this.main_data.filter(d => d.Well==well)[0];
  }

  extract_fluorescence_data(amp_data, well, target) {
    /**
     * Extracts the raw Rn values from the amplification data
     * for one well and one target. This is a very inneficient way to
     * do this, but I think it will still be fast enough.
     */
    let tmp_data = amp_data.filter(row => ((row.Well==well) && (row['Target Name']==target)));
    let i = 0;
    let fluor_curve = [];
    let cycles_bad = false
    for (let row of tmp_data) {
      i += 1;
      if (row['Cycle']!=String(i)) {
        cycles_bad = true;
      } else {
        fluor_curve.push(row['Delta Rn'].replace(/,/g , "")); // Have to do the replace because THESE MANIACS export data with commas for thousands etc.
      }
    }
    if (cycles_bad) {
      this.plate_errors.push('ERROR: Cycles are out of order in amplification data for well ' + well + ', ' + target);
      return null;
    }
    if (fluor_curve.length!=45) {
      this.plate_warnings.push('WARNING: Unexpected number of cycles (' + String(fluor_curve.length) + ') in amplication data for well' + well + ', ' + target);
    }
    return fluor_curve.join(';')
  }

  make_main_data() {
    /**
     * This function takes the raw_data object from load_file() and creates the main_data dataFrame
     * with these columns:
     * Directly from the results data: "Well", "Well Position"
     * Mapped from the sample setup data: "Sample Name"
     * Added from the results data: "N1 CT", "RDRP CT", "RNASEP CT"
     * Added from the amplification data: "N1 FC", "RDRP FC", "RNASEP FC" 
     * (a semicolon separated list of fluorescence data (FC for flourescence curve))
     * To be added later: 
     * "RawCall" - computed based on CT results
     * "Override" - boolean marks if this has been overriden
     * "OverrideCall" - override call, chosen by user or set by plate failure
     * "FinalCall" - possibly changed by plate failures or overrides
     */
    if (whoo_DEBUG) console.log('in make_main_data()');
    let self = this;
    this.main_data = [];
    this.plate_errors = [];
    this.plate_warnings = [];
    this.CT_cutoffs = {};
    this.sample_name_map = {};
    for (let row of this.raw_data['Sample Setup']) {
      // Getting a map between "Well" (a # identifying a well) and sample names
      // This is a bit redundant because there are three rows per well/sample, but that's ok
      this.sample_name_map[row['Well']] = row['Sample Name']; 
    }
    for (let group of d3.groups(this.raw_data['Results'], d => d.Well)) { // groups data by well number
      let rows = group[1];
      let tmp_obj = {'Well': group[0], 'Well Position': rows[0]['Well Position']};
      for (let row of rows) {
        tmp_obj[row['Target Name'] + ' CT'] = parseFloat(row['CT']); // adding CT values for the appropriate target
        // Adding fluorescence data
        tmp_obj[row['Target Name'] + ' FC'] = this.extract_fluorescence_data(this.raw_data['Amplification Data'], group[0], row['Target Name']); 
        tmp_obj['Sample Name'] = this.sample_name_map[row['Well']];
        // Recording the CT cutoffs  - This is very redundant/inefficient, there are only three key: value pairs
        // but we just keep resetting them as we go through the data.
        // Have to do the replace because THESE MANIACS export data with commas for thousands etc.
        this.CT_cutoffs[row['Target Name']] = Number(row['Ct Threshold'].replace(/,/g , "")); 
      }
      this.main_data.push(tmp_obj)
    }
    this.check_plate_data(this.raw_data['Sample Setup']);
    this.call_wells();
    this.filtered_data = this.main_data;
    this.update_filtered_wells()
    this.make_table();
    this.make_summary();
    this.make_svg();
    d3.select("#whoo_export").on('click', function() { self.export(); });
  }

  check_plate_data(sample_data) {
    /**
     * Checks the control data for the plate
     * - Checks if the control samples are in the expected wells 
     * - Checks if EITHER of both the PCTs failed (were negative, CT>30) for either N1 or RDRP
     * - Checks if BOTH of the NTCs failed (were positive, CT<38) for either N1 or RDRP
     * All of these checks are recorded in an array (this.plate_failures)
     */
    if (whoo_DEBUG) console.log('in check_plate_data()');
    this.plate_failures = [];
    let ntc_sample_names = ['NTC-W', 'NTC-P'];
    this.ntc_data = this.main_data.filter(row => (ntc_sample_names.indexOf(row['Sample Name']) > -1));
    this.ptc_data = this.main_data.filter(row => (row['Sample Name'] == 'PTC'));
    this.ntc_well_positions = d3.map(this.ntc_data, row => row['Well Position']);
    this.ptc_well_positions = d3.map(this.ptc_data, row => row['Well Position']);
    
    // Checking they are in the expected spots, will raise warnings if not
    if (!((this.ptc_well_positions.indexOf('O23')>-1) && (this.ptc_well_positions.indexOf('P24')>-1) && (this.ptc_well_positions.length==2))) {
      this.plate_warnings.push('WARNING: Unexpected plate layout, PTCs are at Wells: ' + this.ptc_well_positions.join(', '));
    }
    if (!((this.ntc_well_positions.indexOf('O24')>-1) && (this.ntc_well_positions.indexOf('P23')>-1) && (this.ntc_well_positions.length==2))) {
      this.plate_warnings.push('WARNING: Unexpected plate layout, NTCs are at Wells: ' + this.ntc_well_positions.join(', '));
    }
    
    // Checking for plate errors
    if (!((this.ptc_well_positions.length==2) && (this.ntc_well_positions.length==2))) {
      this.plate_errors.push('PLATE ERROR: Could not find 2 PTCs and/or 2 NTCs'); // Did not find control wells...
    } else {
      // Note: "Undetermined" CTs are NaN here. NaN compared to a number will always return false
      // Therefore all comparisons with CTs will ask if a CT is less than some threshold
      // CHECK 1: PTCs are NOT both NEG for N1
      if ( (!(this.ptc_data[0]['N1 CT']<30)) && (!(this.ptc_data[1]['N1 CT']<30)) ) this.plate_errors.push('PLATE ERROR: Both PTCs N1 CT>30');
      // CHECK 2: PTCs are NOT both NEG for RDRP
      if ( (!(this.ptc_data[0]['RDRP CT']<30)) && (!(this.ptc_data[1]['RDRP CT']<30)) ) this.plate_errors.push('PLATE ERROR: Both PTCs RDRP CT>30');
      // CHECK 3: NTCs are both NOT POS for N1
      if ((this.ntc_data[0]['N1 CT']<38) || (this.ntc_data[1]['N1 CT']<38)) this.plate_errors.push('PLATE ERROR: At least one NTC N1 CT<38'); 
      // CHECK 3: NTCs are both NOT POS for RDRP
      if ((this.ntc_data[0]['RDRP CT']<38) || (this.ntc_data[1]['RDRP CT']<38)) this.plate_errors.push('PLATE ERROR: At least one NTC RDRP CT<38'); 
    }
  }

  call_wells() {
    /**
     * Calls each well as Positive, Negative, Inconclusive, Rerun-A, or Rerun-B
     */
    if (whoo_DEBUG) console.log('in call_wells()');
    console.log(this.plate_failures, this.loading_errors);
    for (let d of this.main_data) {
      // Note: "Undetermined" CTs are NaN here. NaN compared to a number will always return false
      // Therefore all comparisons with CTs will ask if a CT is less than some threshold
      if ( (d['N1 CT']<30) || (d['RDRP CT']<30) ) { 
        d['RawCall'] = 'Positive'; // If N1 CT < 30 OR RDRP CT < 30, call Positive
      } else if ( (d['N1 CT']<38) && (d['RDRP CT']<38) ) {
        d['RawCall'] = 'Positive'; // If N1 CT < 38 AND RDRP CT < 38, call Positive
      } else if ( ( (!(d['N1 CT']<38)) && (!(d['RDRP CT']<38)) ) && (d['RNASEP CT']<34) ) {
        d['RawCall'] = 'Negative'; // If N1 CT > 38 AND RDRP CT > 38 AND RNASEP CT < 34, call Negative
      } else if ( (d['N1 CT']<38) || (d['RDRP CT']<38) ) {
        d['RawCall'] = 'Rerun-B'; // If N1 CT < 38 OR RDRP CT < 38, call Rerun-B
      } else {
        d['RawCall'] = 'Rerun-A'; // ELSE call Rerun-A
      }
      if (this.plate_errors.length>0) {     // If there are any plate-failure errors, set everything to Rerun-A (by override)
        d['Override'] = true;
        d['OverrideCall'] = 'Rerun-A';
        d['FinalCall'] = 'Rerun-A';    // FinalCall will be equal to RawCall when Override is false, and OverrideCall when Override is true
      } else {                              // ELSE:
        d['Override'] = false;
        d['OverrideCall'] = d['RawCall'];   // Override call defaults to the RawCall, can be changed by the user if Override is true
        d['FinalCall'] = d['RawCall'];      // FinalCall will be equal to RawCall when Override is false, and OverrideCall when Override is true
      }
      d['whoo_class_base'] = this.get_whoo_class_base(d);
    }
    if (whoo_DEBUG) console.log('Data after calls:', d3.groups(this.main_data, d => d.RawCall));
  }

  //#endregion

  get_whoo_class_base(d) {
    /**
     * Given a data row, returns the current styling class additions for that row:
     * whoo_control: only if it is a control
     * raw_call_CALL: a class for the raw call (e.g. raw_call_Positive)
     * ...
     * there will be other class styling, but these are the immutable parts
     */
    let class_str = 'raw_call_'+d['RawCall'];
    if ((this.ntc_well_positions.indexOf(d['Well Position']) > -1) || (this.ptc_well_positions.indexOf(d['Well Position']) > -1)) class_str += ' whoo_control';
    if (this.ntc_well_positions.indexOf(d['Well Position']) > -1) class_str += ' whoo_ntc';
    if (this.ptc_well_positions.indexOf(d['Well Position']) > -1) class_str += ' whoo_ptc';
    return class_str;
  }

  make_table() {
    /**
     * Constructs the table from the main data
     */
    let self = this;
    d3.select('#whoo_table_body').selectAll('.tr')
      .data(this.main_data)
      .enter()
      .append('tr')
        .attr('class', function(d) { return 'whoo_item whoo_table_row ' + d['whoo_class_base']; })
        .style('background-color', d => whoo_colors[d.FinalCall])
        .on('mouseover', function(event, d) { d3.selectAll('.svg_data').classed('hovered_data', td => td.Well==d.Well); }) //hover on table -> hover display on svg
        .on('mouseout', function(event, d) { d3.selectAll('.svg_data').classed('hovered_data', false); }) // nothing hovered
        .on('click', function(event, d) { self.highlight_well(d.Well); });

    d3.selectAll('.whoo_table_row').append('td').attr('class', 'Well_row_entry').html(d => d['Well Position']);
    d3.selectAll('.whoo_table_row').append('td').attr('class', 'SampleName_row_entry').html(d => d['Sample Name']);
    d3.selectAll('.whoo_table_row').append('td').attr('class', 'RawCall_row_entry').html(d => d.RawCall);
    let override_entry = d3.selectAll('.whoo_table_row').append('td').attr('class', 'Override_row_entry');
    d3.selectAll('.whoo_table_row').append('td').attr('class', 'CT_row_entry').html(d => (d['N1 CT']) ? d['N1 CT'].toFixed(1) : 'Und.'); // Or (||) override makes NaN show up as -
    d3.selectAll('.whoo_table_row').append('td').attr('class', 'CT_row_entry').html(d => (d['RDRP CT']) ? d['RDRP CT'].toFixed(1) : 'Und.');
    d3.selectAll('.whoo_table_row').append('td').attr('class', 'CT_row_entry').html(d => (d['RNASEP CT']) ? d['RNASEP CT'].toFixed(1) : 'Und.');
    override_entry.append('input')
      .attr('type', 'checkbox')
      .attr('class', 'override_checkbox')
      .property('checked', d => d.Override)
      .on('change', function(event, d) { self.override_change(d, this.parentNode); }); // checkbox change triggers function call
    let override_select = override_entry.append('select')
      .attr('class', 'override_select')
      .style('display', d => (d.Override) ? 'inline' : 'none')
      .on('change', function(event, d) {  
        d.OverrideCall = this.value;     // select box change triggers OverrideCall to change
        d.FinalCall = d.OverrideCall;    // and FinalCall to change to match it
        d3.select(this.parentNode.parentNode).style('background-color', d => whoo_colors[d.FinalCall]);
        self.update_calls_display();
      });
    for (let res of whoo_possible_results) {
      override_select.append('option').attr('value', res).html(res);
    }
    override_select.property('value', d => d.OverrideCall);

    d3.select('#unoverride').on('click', function() { self.set_override_all(false); }); // sets up listener for unoverride button press
    d3.select('#fail_plate').on('click', function() { self.set_override_all(true); }); // sets up listener for fail plate button press

    d3.select('#prev_well').on('click', function() { self.iterate_over_wells(-1); }); // sets up listener for previous button press
    d3.select('#next_well').on('click', function() { self.iterate_over_wells(1); });  // sets up listener for next button press
  }

  update_filtered_wells() {
    /**
     * updates this.filtered_well_nums, which holds an array of all wells in the currently filtered set
     */
    this.filtered_well_nums = [];
    for (let d of this.filtered_data) {
      this.filtered_well_nums.push(Number(d.Well));
    }
  }

  filter_by_call(call) {
    /**
     * On click on one of the calls in the summary, filter the table
     */
    if (whoo_DEBUG) console.log('in filter_by_call()', call);
    let self = this;
    if (call == 'show_all') {
      d3.selectAll('.whoo_table_row').classed('filtered_out_data', false);
      self.filtered_data = self.main_data;
    } else {
      // Deprecated: if you clicked on the already highlighted well, unhighlight, otherwise highlight the well you clicked on
      // self.filtered_call = (self.filtered_call==call) ? null : call;  DECIDED NOT TO DO THIS BEHAVIOR SINCE I HAVE A SHOW ALL BUTTON
      self.filtered_call = call;
      if (self.filtered_call) {
        if (self.filtered_call == 'controls') {
          d3.selectAll('.whoo_table_row').classed('filtered_out_data', d => control_names.indexOf(d['Sample Name']) == -1);
          self.filtered_data = self.main_data.filter(d => control_names.indexOf(d['Sample Name']) > -1);
        } else {
          d3.selectAll('.whoo_table_row').classed('filtered_out_data', d => ((d.FinalCall!=self.filtered_call) || (control_names.indexOf(d['Sample Name']) > -1)));
          self.filtered_data = self.main_data.filter(d => ((d.FinalCall==self.filtered_call) && (control_names.indexOf(d['Sample Name']) == -1)));
        }
      } else {
        d3.selectAll('.whoo_table_row').classed('filtered_out_data', false);
        self.filtered_data = self.main_data;
      }
    }
    this.update_filtered_wells()
  }

  make_summary() {
    /**
     * Adds buttons for each result type to the summary panel
     * determines how many are in each FinalCall category by filtering the main dataset (controls excluded)
     */
    let self = this;
    d3.select('#summary_display').selectAll('.result_button')
      .data(whoo_possible_results)
      .enter()
      .append('div')
        .attr('class', 'whoo_button result_button summary_result')
        .style('background-color', d => whoo_colors[d])
        .html(function(d) { return String(self.main_data.filter(td => ((td.FinalCall==d) && (control_names.indexOf(td['Sample Name']) == -1))).length) + ' ' + d; })
        .on('click', function(event, d) { self.filter_by_call(d); });
    
    d3.select('#summary_display')
      .append('div')
        .attr('class', 'whoo_button result_button')
        .html('Controls')
        .on('click', function() { self.filter_by_call('controls'); });
    
    d3.select('#summary_display')
      .append('div')
        .attr('class', 'whoo_button result_button')
        .html('Show All')
        .on('click', function() { self.filter_by_call('show_all'); });
  }

  iterate_over_wells(step) {
    /**
     * highlights the previous (step=-1) or next (step=1) well in the set of filtered wells
     */
    if (this.filtered_well_nums.length>0) {
      let current_well = Number(this.highlighted_well) || 0; // or operator override sets it to 0 if it is null
      current_well += step;
      while (this.filtered_well_nums.indexOf(current_well)==-1) {
        if (current_well < 0) {
          current_well = 384;
        } else if (current_well > 384) {
          current_well = 1;
        } else {
          current_well += step
        }
      }
      this.highlight_well(String(current_well));
    }
  }

  scroll_to_row(row_element) {
    /**
     * Checks if the highlighted well row is in the display, if not scrolls so it is at the top
     */
    let table_holder = d3.select('#whoo_table_holder');
    let div_height = table_holder.node().getBoundingClientRect().height;
    if (row_element.offsetTop >= div_height+table_holder.property('scrollTop')) {
      table_holder.property('scrollTop', row_element.offsetTop-whoo_table_row_height);
    } else if (row_element.offsetTop < table_holder.property('scrollTop')+whoo_table_row_height) {
      table_holder.property('scrollTop', row_element.offsetTop-div_height+whoo_table_row_height+5);
    }
  }

  highlight_well(well) {
    /**
     * On click on the svg or the table, highlights on the svg and the table
     */
    // if you clicked on the already highlighted well, unhighlight, otherwise highlight the well you clicked on
    if (whoo_DEBUG) console.log('in highlight_well()', well);
    let self = this;
    self.highlighted_well = (self.highlighted_well==well) ? null : well; 
    d3.selectAll('.svg_data').classed('clicked_data', d => d.Well==self.highlighted_well);
    let row_element;
    let well_position = '';
    d3.selectAll('.whoo_table_row').classed('clicked_data', function(d) {
      if (d.Well==self.highlighted_well) {
        row_element = this; // using this classed call to also find the row element
        well_position = d['Well Position']; // And grab the well position
      }
      return (d.Well==self.highlighted_well);
    });
    d3.select('#well_display').html(well_position);
    this.scroll_to_row(row_element);
  }

  make_svg() {
    /**
     * Makes the data browser svg
     */
    let self = this;
    this.highlighted_well = null;
    this.svg = d3.select('#whoo_svg').append('g').attr('class', 'whoo_item'); // one g to hold them all
    this.cycle_scale = d3.scaleLinear().domain([0,45]).range([100,550]);
    this.fluor_scale = d3.scaleLog().domain([100,1000000]).range([285,50]);
    this.axis_x_scale = d3.scaleLinear().domain([0,1]).range([100,550]);
    this.axis_y_scale = d3.scaleLinear().domain([0,1]).range([285,50]);
    this.fline = d3.line()
      .x(function(d) { return self.cycle_scale(d.x); })
      .y(function(d) { return (d.y > 100) ? self.fluor_scale(d.y) : self.fluor_scale(100); }); // clipping numbers <1 for log scale

    this.axis_line = d3.line()
      .x(function(d) { return self.axis_x_scale(d.x); })
      .y(function(d) { return self.axis_y_scale(d.y); });

    // Making cutoff lines
    let i = 0;
    for (let gene of Object.keys(whoo_gene_colors)) {
      this.svg.append('path')
        .attr('class', 'ct_cutoff')
        .attr('stroke', whoo_gene_colors[gene])
        .attr('d', function(d) {
          return self.fline([{'x': 1, 'y': self.CT_cutoffs[gene]}, {'x': 45, 'y': self.CT_cutoffs[gene]}]);
        });

      this.svg.append('path')
        .attr('class', 'legend_line')
        .attr('stroke', whoo_gene_colors[gene])
        .attr('d', function(d) {
          return self.axis_line([{'x': i*0.3+0.1, 'y': 1.1}, {'x': i*0.3+0.16, 'y': 1.1}]);
        });

      this.svg.append('text')
        .attr('class', 'legend_text')
        .style('fill', whoo_gene_colors[gene])
        .style('text-align', 'start')
        .style('alignment-baseline', 'middle')
        .attr('x', self.axis_x_scale(i*0.3+0.18))
        .attr('y', self.axis_y_scale(1.1))
        .html(gene);
      i += 1;
    }

    let row_range = [50, 285]; // defining a variable in order to offset scaleBand locations so ticks match...
    this.column_scale = d3.scaleLinear().domain([1,24]).range([650,1100]);
    this.row_scale = d3.scaleBand().domain(whoo_rows).range(row_range);

    this.svg.selectAll('.svg_data')
      .data(self.main_data)
      .enter()
      .append('g')
        .attr('class', function(d) { return 'svg_data ' + d['whoo_class_base']; })
        .on('click', function(event, d) { self.highlight_well(d.Well); });

    for (let gene of Object.keys(whoo_gene_colors)) {
      this.svg.selectAll('.svg_data')
        .append('path')
        .attr('class', 'fluor_path')
        .attr('stroke', whoo_gene_colors[gene])
        .attr('d', function(d) {
          let yvals = d[gene + ' FC'].split(';');
          let tmp_line_data = [];
          for (let i=0; i<yvals.length; i++) {
            tmp_line_data.push({'y': Number(yvals[i]), 'x': i+1});
          }
          return self.fline(tmp_line_data);
        });
    }
    
    this.svg.append('g').attr("transform", "translate(0, "+self.fluor_scale(100)+")").call(d3.axisBottom().scale(self.cycle_scale));
    this.svg.append('g').attr("transform", "translate("+String(self.cycle_scale(0))+", 0)").call(d3.axisLeft().scale(self.fluor_scale));
    this.svg.append('text').html('Cycle')
      .attr('text-anchor', 'middle')
      .attr('class', 'axis_label')
      .attr('x', self.cycle_scale(45/2))
      .attr('y', self.fluor_scale(100)+40);

    this.svg.append('text').html('\u0394 Rn')
      .attr('text-anchor', 'end')
      .attr('class', 'axis_label')
      .attr('x', self.cycle_scale(0)-40)
      .attr('y', self.fluor_scale(20000));

    this.svg.append('g')
      .attr('class', 'plate_axis')
      .attr("transform", "translate("+String(self.column_scale(0.5))+", 0)")
      .call(d3.axisLeft().scale(self.row_scale));
    this.svg.append('g')
      .attr('class', 'plate_axis')
      .attr("transform", "translate(0, "+String(self.row_scale('A'))+")")
      .call(d3.axisTop().scale(self.column_scale).ticks(24));

    this.svg.selectAll('.svg_data')
      .append('circle')
      .attr('class', 'well_icon')
        .attr('cx', d => self.column_scale(parseInt(d['Well Position'].slice(1,))))
        .attr('cy', d => self.row_scale(d['Well Position'][0])+(row_range[1]-row_range[0])/32)
        .attr('r', 6)
        .attr('stroke', '#EEE')
        .attr('fill', d => whoo_colors[d.FinalCall].replace('33%', '45%'));

    d3.select('#show_controls')
      .on('click', function() {
        let controls_button = d3.select(this);
        if (controls_button.classed('show_controls_on')) {
          controls_button.classed('show_controls_on', false);
          d3.selectAll('.whoo_control').classed('whoo_control_on', false);
          controls_button.html('show controls');
        } else {
          controls_button.classed('show_controls_on', true);
          d3.selectAll('.whoo_control').classed('whoo_control_on', true);
          controls_button.html('hide controls');
        }
      })

  }

  update_calls_display() {
    /**
     * updates all display elements after any call has changed
     * updates summary counts
     * determines how many are in each FinalCall category by filtering the main dataset (controls excluded)
     * updates colors of well icons on the svg
     */
    let self = this;
    d3.select('#summary_display').selectAll('.summary_result')
      .html(function(d) { return String(self.main_data.filter(td => ((td.FinalCall==d) && (control_names.indexOf(td['Sample Name']) == -1))).length) + ' ' + d; });

    this.svg.selectAll('.well_icon')
        .attr('fill', d => whoo_colors[d.FinalCall].replace('33%', '45%') );
  }

  set_override_all(override_bool) {
    /**
     * Sets the whole plate override settings
     * false:
     * Reverses all overrides - this may be used when a whole plate is declared inconclusive
     * due to a problem with controls, but the user decides the controls are ok and wants to
     * use the raw calls
     * true:
     * Sets all wells to by overridden and Inconclusive
     */
    if (whoo_DEBUG) console.log('in set_override_all()');
    d3.selectAll('.override_checkbox')
      .each( function(d) {
        d.Override = override_bool;        // set overrides true or false
        d.OverrideCall = (override_bool) ? 'Inconclusive' : d.RawCall; // initialize override calls to Inconclusive or RawCall (true or false respectively)
        d.FinalCall = (override_bool) ? d.OverrideCall : d.RawCall;   // set FinalCall to OverrideCall or RawCall (true or false respectively)
        d3.select(this.parentNode).selectAll('.override_select')
          .style('display', (override_bool) ? 'inline' : 'none') // show or hide select
          .property('value', d.OverrideCall); // set default value to OverrideCall value
        d3.select(this.parentNode.parentNode).style('background-color', td => whoo_colors[td.FinalCall]); // set color to new FinalCall
        d3.select(this).property('checked', override_bool); // display override check correctly
      });
    this.update_calls_display();
  }

  override_change(d, row_element) {
    /**
     * This function gets called when an override checkbox is clicked
     * If the override is already on, this:
     *  - sets the Override property to false
     *  - hides the corresponding override_checkbox
     *  - sets the FinalCall property to RawCall
     */
    if (whoo_DEBUG) console.log('in override_change()', d);
    if (d.Override) {
      d.Override = false;
      d.FinalCall = d.RawCall;
      d3.select(row_element).selectAll('.override_select').style('display', 'none');
    } else {
      d.Override = true;
      d.FinalCall = d.OverrideCall;
      d3.select(row_element).selectAll('.override_select').style('display', 'inline');
    }
    d3.select(row_element.parentNode).style('background-color', d => whoo_colors[d.FinalCall]);
    this.update_calls_display();
  }
}

function setup() {
  //parse_QS7F_txt_output('test_data/2021-02-17 135556.txt');
}

function new_file_handler() {
  d3.select("#all_contents").style('display', 'block');
  whoo_holder = new whooData();
  whoo_holder.load_file();
}