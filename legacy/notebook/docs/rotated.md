---
title: Rotating genomes
---

# Alignment lengths

HBV is a circular genome and so the genomes found on NCBI have been linearised.
The following plot shows the % of the genome which aligns to the reference (in black), where around 1000 genomes seem to only have 50% matching.
If we rotate the genomes in order to all start at the same place then we can improve this (coloured dots).

Genomes with over 90% match: ${numAbove(rotatedMatchCount, 90)} (rotated) vs ${numAbove(rawMatchCount, 90)} (NCBI)

Genomes with over 80% match: ${numAbove(rotatedMatchCount, 80)} vs ${numAbove(rawMatchCount, 80)}

```js
display(Plot.plot({
  marginTop: 20,
  marginRight: 20,
  marginBottom: 30,
  marginLeft: 60,
  grid: true,
  inset: 10,
  // aspectRatio: 0.6,
  width: 1000,
  figure: true,
  color: {legend: true},
  className: "largerFont",
  x: { tickSize: 15, label: "Genomes"},
  y: { tickSize: 15, label: "Aln match (%)"},
  marks: [
    Plot.frame(),
    Plot.dot(rawMatchCount, {x: "idx", y: "count", stroke: "black"}),
    Plot.dot(rotatedMatchCount, {x: "idx", y: "count", stroke: "genotype", fill: "genotype", symbol: "circle", r: 3}),
    Plot.crosshair(rotatedMatchCount, {x: "idx", y: "count", stroke: "genotype"}),
  ]
}))
```
<style type="text/css">
  .largerFont {
    font-size: 16px;
  }
</style>



```js
const parseMetadata = async () => {
    return Object.fromEntries(
        (await FileAttachment("data/metadata.tsv").tsv())
            .map((row) => [row.name, row])
    )
}
const metadata = await parseMetadata()
```



```js
const refLength = 3182;
const parseAlignment = async (attachment) => {
  return Object.entries(
      await attachment.json()
    )
    .map(([name, count]) => ({name, count}))
    .sort((a, b) => a.count > b.count ? -1 : 1)
    .map((d, idx) => ({...d, count: d.count/refLength*100}))
    .map((d, idx) => ({...d, idx}))
    .map((d) => ({...d, genotype: metadata[d.name]?.genotype_genbank}))
}
```

```js
const rotatedMatchCount = await parseAlignment(FileAttachment("data/alignment_rotated.json"))
const rawMatchCount = await parseAlignment(FileAttachment("data/alignment_raw.json"))
// console.log("rotatedMatchCount", rotatedMatchCount.slice(0, 10))
// console.log("rawMatchCount", rawMatchCount.slice(0, 10))
```



```js
const numAbove = (dataset, perc) => {
  return dataset.filter((d) => d.count > perc).length
}
```