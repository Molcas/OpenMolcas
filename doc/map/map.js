/* helper function to run code when an image is loaded */
async function loadImage(imageUrl) {
  let img;
  const imageLoadPromise = new Promise(resolve => {
    img = new Image();
    img.onload = resolve;
    img.src = imageUrl;
  });
  await imageLoadPromise;
  return img;
}

/* hackish function to add a data-url shadow image based on the icon */
function add_shadow(group, icon) {
  var opacity = 0.3;
  var angle = 35;
  var scale = 0.5;
  var blur = 1.5;
  var blurscale = 2*blur;
  var tan = Math.tan(angle*Math.PI/180);
  var canvas = document.createElement('canvas');
  var context = canvas.getContext('2d');
  loadImage(icon.options.iconUrl).then(img =>  {
    var w = img.naturalWidth;
    var h = img.naturalHeight;
    canvas.width = w+h*tan+2*blurscale;
    canvas.height = h*scale+2*blurscale;
    context.setTransform(1,0,-tan,scale,0,0);
    context.filter = 'blur('+blur+'px)';
    context.drawImage(img, (h+blurscale)*tan+blurscale, blurscale/scale);
    context.globalCompositeOperation = 'source-in';
    context.fillStyle = 'rgba(0,0,0,'+opacity+')';
    context.fillRect(0,0,(w+h*tan+blurscale)/scale,h/scale);
    var pngUrl = canvas.toDataURL('image/png');
    icon.options.shadowUrl = pngUrl;
    icon.options.shadowSize = [canvas.width, canvas.height];
    icon.options.shadowAnchor = [devIcon.options.iconAnchor[0]+blurscale, canvas.height-blurscale];
    group.eachLayer(function(layer) {
      layer.setIcon(icon);
    });
  });
}

/* ditto to colorize icon */
function color_icon(group, icon, color) {
  var canvas = document.createElement('canvas');
  var context = canvas.getContext('2d');
  loadImage(icon.options.iconUrl).then(img =>  {
    var w = img.naturalWidth;
    var h = img.naturalHeight;
    canvas.width = w;
    canvas.height = h;
    context.drawImage(img,0,0);
    context.globalCompositeOperation = 'source-in';
    context.fillStyle = color;
    context.fillRect(0,0,w,h);
    var pngUrl = canvas.toDataURL('image/png');
    icon.options.iconUrl = pngUrl;
    group.eachLayer(function(layer) {
      layer.setIcon(icon);
    });
  });
}

/* function to create a layer from an array of markers */
var markers = L.featureGroup();
function markerLayer(inMarkers, icon) {
  var outMarkers = L.layerGroup();
  for (var i = 0; i < inMarkers.length; ++i) {
    var m = L.marker([inMarkers[i].lat, inMarkers[i].lng], {icon: icon})
      .bindPopup('<a href="' + inMarkers[i].url + '" target="_blank">' + inMarkers[i].name + '</a><p>' + inMarkers[i].city + '</p>');
    if (inMarkers[i].virt) { m.setOpacity(0.5); }
    outMarkers.addLayer(m);
    markers.addLayer(m);
    // add periodic copies
    var m = L.marker([inMarkers[i].lat, inMarkers[i].lng+360], {icon: icon})
        .bindPopup('<a href="' + inMarkers[i].url + '" target="_blank">' + inMarkers[i].name + '</a><p>' + inMarkers[i].city + '</p>');
    if (inMarkers[i].virt) { m.setOpacity(0.5); }
    outMarkers.addLayer(m);
    var m = L.marker([inMarkers[i].lat, inMarkers[i].lng-360], {icon: icon})
      .bindPopup('<a href="' + inMarkers[i].url + '" target="_blank">' + inMarkers[i].name + '</a><p>' + inMarkers[i].city + '</p>');
    if (inMarkers[i].virt) { m.setOpacity(0.5); }
    outMarkers.addLayer(m);
  }
  add_shadow(outMarkers, icon);
  return outMarkers;
};

var molcasIcon = L.Icon.extend({
  options: {
    iconSize: [50,46],
    iconAnchor: [18,46],
    popupAnchor: [0,-46],
  }
});
var devIcon = new molcasIcon({iconUrl: 'openmolcas.png'})
var workIcon = new molcasIcon({iconUrl: '50molcas.png'})

var osm = L.tileLayer('http://tile.openstreetmap.org/{z}/{x}/{y}.png', {
  attribution: '© <a href="http://www.openstreetmap.org/copyright" target="_blank">OpenStreetMap</a>',
  maxZoom: 19,
});

var otm = L.tileLayer('https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png', {
  attribution: 'Map data: © <a href="https://openstreetmap.org/copyright" target="_blank">OpenStreetMap</a> contributors, <a href="viewfinderpanoramas.org" target="_blank">SRTM</a> | Map style: © <a href="https://opentopomap.org" target="_blank">OpenTopoMap</a> (<a href="https://creativecommons.org/licenses/by-sa/3.0/" target="_blank">CC-BY-SA</a>)',
  maxZoom: 17,
});

var workshops = markerLayer(dev_workshops, workIcon);
color_icon(workshops, workIcon, '#010165');

var institutions = markerLayer(dev_places, devIcon);

var map = L.map('map', {
  center: [25.0, 5.0],
  zoomSnap: 0.05,
  minZoom: 2,
  zoom: 2.75,
  layers: [osm, institutions],
  worldCopyJump: true,
});

L.control.scale().addTo(map);

var baseMaps = {
  'OpenStreetMap': osm,
  'OpenTopoMap': otm,
};

var overlayMaps = {
  'Workshops': workshops,
  'Institutions': institutions,
};

var layerControl = L.control.layers(baseMaps, overlayMaps).addTo(map);

map.fitBounds(markers.getBounds(), {padding: [50, 50]});
