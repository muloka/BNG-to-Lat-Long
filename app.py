import os
import json
from flask import Flask
from LatLongUTMConversion import bng_to_latlon

app = Flask(__name__)


@app.route('/bng/')
@app.route('/')
def index():
    return 'To use API visit http://[host]/bng/[easting]/[northing]'


@app.route('/bng/<int:easting>/<int:northing>')
def bng(easting, northing):
    coordinates = bng_to_latlon(easting, northing)
    tmp = {"lat": coordinates[0], "long": coordinates[1]}
    return json.dumps(tmp)


if __name__ == '__main__':
    # Bind to PORT if defined, otherwise default to 5000.
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)
