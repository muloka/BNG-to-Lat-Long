import os
import json
from flask import Flask
from osgeo.osr import SpatialReference, CoordinateTransformation

# Define BDA2000 / Bermuda 2000 National Grid (EPSG 3770)
epsg3770 = SpatialReference()
epsg3770.ImportFromEPSG(3770)

# Define the wgs84 system (EPSG 4326)
epsg4326 = SpatialReference()
epsg4326.ImportFromEPSG(4326)

bng2latlon = CoordinateTransformation(epsg3770, epsg4326)

app = Flask(__name__)

@app.route('/bng/')
@app.route('/')
def index():
    return 'To use API visit http://[host]/bng/[easting]/[northing]'


@app.route('/bng/<int:easting>/<int:northing>')
def bng(easting, northing):
    coordinates = bng2latlon.TransformPoint(easting, northing)
    tmp = { "lat": coordinates[1], "lon": coordinates[0] }
    return json.dumps(tmp)


if __name__ == '__main__':
    # Bind to PORT if defined, otherwise default to 5000.
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)
