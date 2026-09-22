# Terminal artifacts emitted after sealing

As in v52, the destructive POSTSEAL report, atomic completion token, and delivery index are emitted only after the deterministic bundle has been sealed and replayed. They are therefore distributed next to the bundle rather than embedded in it.
