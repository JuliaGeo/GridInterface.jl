```@meta
CurrentModule = GridInterface
```

# Subsets, chunks and wrappers

## Subsets and chunks

A subset is a grid that remembers its parent; a chunking is a declared partition into subsets, each with a cheap covering extent.

```@docs
Subset
subset
parentpositions
union
grow
chunks
chunkat
LevelChunks
hierarchy_chunks
AncestorChunks
WindowChunks
```

## Wrappers

A wrapper is a grid plus something the caller holds for it. Positions are unchanged, so every verb forwards to the parent and the wrapper overrides the one or two it changes. GridInterface holds no global state: these wrappers are its only caches.

```@docs
AbstractWrappedGrid
WRAPPED_VERBS
GridAndTree
ChunkedGrid
mockchunks
GridAndVertices
```
