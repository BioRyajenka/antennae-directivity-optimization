class NDArrayEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            output = io.BytesIO()
            np.savez_compressed(output, obj=obj)
            return {'b64npz' : base64.b64encode(output.getvalue())}
        return json.JSONEncoder.default(self, obj)


def ndarray_decoder(dct):
    if isinstance(dct, dict) and 'b64npz' in dct:
        output = io.BytesIO(base64.b64decode(dct['b64npz']))
        output.seek(0)
        return np.load(output)['obj']
    return dct

# Make expected non-contiguous structured array:
expected = np.arange(10)[::2]
expected = expected.view('<i4,<f4')

dumped = json.dumps(expected, cls=NDArrayEncoder)
result = json.loads(dumped, object_hook=ndarray_decoder)

assert result.dtype == expected.dtype, "Wrong Type"
assert result.shape == expected.shape, "Wrong Shape"
assert np.array_equal(expected, result), "Wrong Values"