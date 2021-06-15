let point = (x, y, z) => -x * 1e032 - y * 1e013 - z * 1e021 + 1e123;
let ssqrt = (x) => (1 + x).Normalized
let sqrtn = (x) => ((1 + x) * (1 + x.Grade(1) - 0.5 * x.Grade(4))).Normalized;
let sqrt = (x) => (1 + x) / (Math.sqrt(2 + 2 * x[0])) * (1 - x.Grade(4) / (2 + 2 * x.Grade(0)))
function normalize(x, y, z) {
	let norm = Math.sqrt(x * x + y * y + z * z);
	return [x / norm, y / norm, z / norm];
}
function plane(d, n1_, n2_, n3_) {
	let [n1, n2, n3] = normalize(n1_, n2_, n3_);
	return d * 1e0 + n1 * 1e1 + n2 * 1e2 + n3 * 1e3;
}