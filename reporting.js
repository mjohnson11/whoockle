
var main_dat;
var use_dat;

var organizations;

var possible_results = ['NEGATIVE', 'INCONCLUSIVE RESULT', 'POSITIVE'];
var pr_full = [];
for (let p of possible_results) {
  pr_full.push(p+' for SARS-CoV-2');
}
var res_colors = {'POSITIVE': 'hsl(350, 100%, 66%)', 'NEGATIVE': 'hsl(0, 0%, 66%)', 'INCONCLUSIVE RESULT': 'hsl(34, 100%, 66%)'};

var x_scale = d3.scaleLinear().domain([0,1]).range([0,350]);

var full_date_range;
var current_date_range;
var day_num_to_date_string = {}; //hacky way to convert back, just storing the strings

var grid;
var columns = [
  {id: "org", name: "Org", field: "Organization name", width: 200},
  {id: "res", name: "Result", field: "Result text", width: 200},
  {id: "date", name: "Date", field: 'Sample date authorised', width: 200},
  {id: "bc", name: "BC", field: 'Sample sample vial id', width: 200}
];
var options = {
  enableCellNavigation: true,
  enableColumnReorder: false
};

var clicked_result = '';

function click_on_result(res) {
  if (res == clicked_result) {
    clicked_result = '';
  } else {
    clicked_result = res;
  }
  d3.selectAll('results_div').classed('highlighted_res_div', d => d==clicked_result);
  let tmp_dat;
  if (clicked_result == '') {
    tmp_dat = use_dat;
  } else {
    tmp_dat = use_dat.filter(d => d['Result text']==clicked_result+' for SARS-CoV-2');
  }
  grid = new Slick.Grid("#results_table", tmp_dat, columns, options);
}

function date_number(date) {
  // input: date like YYYY-MM-DD
  // output: number like YYYY0MM0DD
  let date_fields = date.split('-');
  return parseInt(date_fields[0])*1000000+parseInt(date_fields[1])*1000+parseInt(date_fields[2]);
}

function parse_date(row) {
  // date strings in sample manager are like "4/12/2021 10:00 AM"
  // I will change them to day only dates in standard date string format https://tc39.es/ecma262/#sec-date-time-string-format
  // and then parse them
  let date_string = row['Sample date authorised'];
  let american_date = date_string.split(' ')[0];
  let date_fields = american_date.split('/');
  let date_string_fixed = date_fields[2]+'-'+date_fields[0].padStart(2, "0")+'-'+date_fields[1].padStart(2, "0");
  row['day'] = date_number(date_string_fixed);
  day_num_to_date_string[row['day']] = date_string_fixed; // hacky, see note above
}

function show_the_stats() {
  let orgs_checked = [];
  d3.selectAll('.org_check_input').each(function(d) { 
    if (d3.select(this).property('checked')) {
      orgs_checked.push(d);
    } 
  })
  console.log('Filter by orgs', orgs_checked);
  console.log('Filter by date', current_date_range);
  current_date_range = [date_number(d3.select('#start_date').property('value')), date_number(d3.select('#end_date').property('value'))];
  use_dat = main_dat.filter(d=>((orgs_checked.indexOf(d['Organization name'])>-1) && (d['day']>=current_date_range[0]) && (d['day']<=current_date_range[1])));
  console.log('Filtered data', use_dat);
  let total = use_dat.length;
  x_scale.domain([0, total]);
  d3.selectAll('.results_bar')
    .style('width', d => String(Math.round(x_scale(use_dat.filter(td=>td['Result text']==d+' for SARS-CoV-2').length)))+'px');
  d3.selectAll('.results_text')
    .style('left', d => String(Math.round(x_scale(use_dat.filter(td=>td['Result text']==d+' for SARS-CoV-2').length))+30)+'px')
    .html(d => String(use_dat.filter(td=>td['Result text']==d+' for SARS-CoV-2').length) + " " + d.split(' ')[0]);

  d3.select("#total_text").html('Total Tests: ' + String(use_dat.length));

  grid = new Slick.Grid("#results_table", use_dat, columns, options);
  
}

function setup() {
  for (let row of main_dat) {
    parse_date(row);
  }
  use_dat = main_dat;
  grid = new Slick.Grid("#results_table", use_dat, columns, options);
  current_date_range = d3.extent(main_dat.map(d=>d['day']));
  full_date_range = current_date_range.map(d => day_num_to_date_string[d]);
  console.log(full_date_range);
  d3.select('#start_date').attr('min', full_date_range[0]).attr('max', full_date_range[1]).attr('value', full_date_range[0]);
  d3.select('#end_date').attr('min', full_date_range[0]).attr('max', full_date_range[1]).attr('value', full_date_range[1]);
  d3.select('#all_contents').style('display', 'block');
  d3.select('#stats_zone').selectAll('.results_div')
    .data(possible_results)
    .enter()
    .append('div')
      .attr('class', 'results_div')
      .on('click', function(e, d) { click_on_result(d); })

  d3.selectAll('.results_div')
    .append('div')
      .attr('class', 'results_bar')
      .style('background-color', d => res_colors[d])
  d3.selectAll('.results_div')
    .append('p')
      .attr('class', 'results_text')
      .style('color', d => res_colors[d])
      .html(d => d.split(' ')[0]);

  d3.selectAll()
  let orgs = new Set(main_dat.map(d=>d['Organization name']));
  organizations = Array.from(orgs);
  d3.select('#org_checkboxes').selectAll('.org_check')
    .data(organizations)
    .enter()
    .append('div')
      .attr('class', 'org_check');
  d3.selectAll('.org_check')
    .append('input')
      .attr('type', 'checkbox')
      .attr('class', 'org_check_input')
      .attr('name', d=>d)
      .property('checked', true)
      .on('change', function() { show_the_stats(); });
  d3.selectAll('.org_check')
    .append('label')
      .attr('class', 'org_check_label')
      .attr('for', d=>d)
      .text(d=>d);
  show_the_stats();
}

function new_stats_file() {
  f_in = document.getElementById('inputfile').files[0];
  let reader = new FileReader();
  reader.readAsText(f_in);
  reader.onload = function() {
      let tmp_dat = d3.csvParse(reader.result, d3.autoType);
      console.log(tmp_dat);
      main_dat = tmp_dat.filter(d=>((d['Sample type']=='COLOR') && (d['Status']=='Authorised') && (d['Component name']=='Interpretive Result') && (pr_full.indexOf(d['Result text'])>-1) && (d['Result sent to color']=="Yes")));
      console.log(main_dat);
      setup();
  }
}

