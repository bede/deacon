export const SEQ_RE = /\.(fastq|fq|fasta|fa)(\.gz)?$/i;

const MATE_PATTERNS = [
  /(^|[._])R([12])(?=([._-]|$))/i,
  /(_)([12])$/,
];

function relativePath(file) {
  return file.webkitRelativePath || file.name;
}

function splitPath(path) {
  const idx = path.lastIndexOf("/");
  if (idx === -1) return { dir: "", name: path };
  return { dir: path.slice(0, idx), name: path.slice(idx + 1) };
}

function stripSeqExt(name) {
  return name.replace(SEQ_RE, "");
}

function mateInfo(file) {
  const rel = relativePath(file);
  const { dir, name } = splitPath(rel);
  const stem = stripSeqExt(name);

  for (const pattern of MATE_PATTERNS) {
    const match = pattern.exec(stem);
    if (!match) continue;
    const mate = match[2] === "1" ? "r1" : "r2";
    const keyStem =
      stem.slice(0, match.index + match[1].length) + "{mate}" + stem.slice(match.index + match[0].length);
    return {
      mate,
      key: `${dir}/${keyStem}`.replace(/^\//, ""),
      dir,
      bareDigit: pattern === MATE_PATTERNS[MATE_PATTERNS.length - 1],
    };
  }
  return null;
}

// Group every "<dir>/<prefix>_" onto the set of trailing integers that follow it,
// so we can tell a true _1/_2 mate pair from ONT read chunks numbered _0, _1, _2, _3, …
function numericSiblingSets(seqs) {
  const sets = new Map();
  for (const file of seqs) {
    const { dir, name } = splitPath(relativePath(file));
    const stem = stripSeqExt(name);
    const match = /_(\d+)$/.exec(stem);
    if (!match) continue;
    const prefixKey = `${dir}/${stem.slice(0, match.index + 1)}`.replace(/^\//, "");
    const set = sets.get(prefixKey) || new Set();
    set.add(Number(match[1]));
    sets.set(prefixKey, set);
  }
  return sets;
}

export function pairSequenceFiles(files) {
  const seqs = files.filter((file) => SEQ_RE.test(file.name));
  const numericSiblings = numericSiblingSets(seqs);
  const pairBuckets = new Map();
  const singles = [];

  for (const file of seqs) {
    const info = mateInfo(file);
    if (!info) {
      singles.push(file);
      continue;
    }
    if (info.bareDigit) {
      // Only pair bare _1/_2 when they are the *only* numbered files sharing this
      // prefix; any other number (_3, _4) means these are treated as ONT.
      const numbers = numericSiblings.get(info.key.replace("{mate}", ""));
      const onlyMates = numbers && [...numbers].every((n) => n === 1 || n === 2);
      if (!onlyMates) {
        singles.push(file);
        continue;
      }
    }
    const bucket = pairBuckets.get(info.key) || { r1: [], r2: [] };
    bucket[info.mate].push(file);
    pairBuckets.set(info.key, bucket);
  }

  const groups = singles.map((file) => ({
    kind: "single",
    file,
    label: relativePath(file),
    size: file.size || 0,
  }));

  for (const [key, bucket] of pairBuckets) {
    if (bucket.r1.length !== 1 || bucket.r2.length !== 1) {
      const seen = [...bucket.r1, ...bucket.r2].map(relativePath).join(", ");
      throw new Error(`Ambiguous or incomplete pair for ${key}: ${seen || "no files"}`);
    }
    const file1 = bucket.r1[0];
    const file2 = bucket.r2[0];
    groups.push({
      kind: "paired",
      file1,
      file2,
      label: `${relativePath(file1)} & ${relativePath(file2)}`,
      size: (file1.size || 0) + (file2.size || 0),
    });
  }

  groups.sort((a, b) => a.label.localeCompare(b.label));
  return groups;
}

// Count work groups the way a user sees them: "1 pair", "2 files", "1 pair and 2 files"
export function describeGroups(groups) {
  const pairs = groups.filter((group) => group.kind === "paired").length;
  const files = groups.length - pairs;
  const parts = [];
  if (pairs) parts.push(`${pairs} pair${pairs === 1 ? "" : "s"}`);
  if (files) parts.push(`${files} file${files === 1 ? "" : "s"}`);
  return parts.join(" and ") || "0 files";
}
