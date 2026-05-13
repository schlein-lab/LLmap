#!/usr/bin/env python3
"""
Generate a minimal ONNX model for testing the FoundationEmbedder.

This creates a simple model that simulates a DNA embedding model:
- Input: token_ids [batch, seq_len] (int64), attention_mask [batch, seq_len] (int64)
- Output: embeddings [batch, embedding_dim] (float32)

The model performs:
1. Cast tokens to float
2. Simple linear projection to embedding space
3. Mean pooling across sequence (respecting attention mask)

This is NOT a real Caduceus-Ph model - just a test surrogate to verify
the ONNX Runtime integration works correctly.
"""

import argparse
import numpy as np
import onnx
from onnx import TensorProto, helper, numpy_helper


def create_test_embedder_model(
    vocab_size: int = 9,
    max_seq_len: int = 514,  # 512 + CLS + SEP
    embedding_dim: int = 256,
    hidden_dim: int = 64,
) -> onnx.ModelProto:
    """Create a minimal test model for DNA embedding."""

    # Input tensors
    input_ids = helper.make_tensor_value_info(
        "input_ids", TensorProto.INT64, ["batch", "seq_len"]
    )
    attention_mask = helper.make_tensor_value_info(
        "attention_mask", TensorProto.INT64, ["batch", "seq_len"]
    )

    # Output tensor
    embeddings = helper.make_tensor_value_info(
        "embeddings", TensorProto.FLOAT, ["batch", embedding_dim]
    )

    # Initialize weights with small random values (deterministic seed)
    np.random.seed(42)

    # Token embedding matrix [vocab_size, hidden_dim]
    token_embed_weights = np.random.randn(vocab_size, hidden_dim).astype(np.float32) * 0.02
    token_embed_init = numpy_helper.from_array(token_embed_weights, name="token_embed_weights")

    # Projection matrix [hidden_dim, embedding_dim]
    proj_weights = np.random.randn(hidden_dim, embedding_dim).astype(np.float32) * 0.02
    proj_init = numpy_helper.from_array(proj_weights, name="proj_weights")

    # Small epsilon for numerical stability
    eps_init = numpy_helper.from_array(np.array([1e-6], dtype=np.float32), name="eps")

    # Build the graph
    nodes = []

    # 1. Gather token embeddings: [batch, seq_len] -> [batch, seq_len, hidden_dim]
    nodes.append(helper.make_node(
        "Gather",
        inputs=["token_embed_weights", "input_ids"],
        outputs=["token_embeddings"],
        axis=0,
    ))

    # 2. Cast attention mask to float: [batch, seq_len]
    nodes.append(helper.make_node(
        "Cast",
        inputs=["attention_mask"],
        outputs=["mask_float"],
        to=TensorProto.FLOAT,
    ))

    # 3. Expand mask for broadcasting: [batch, seq_len, 1]
    nodes.append(helper.make_node(
        "Unsqueeze",
        inputs=["mask_float", "unsqueeze_axes"],
        outputs=["mask_expanded"],
    ))

    # 4. Multiply embeddings by mask
    nodes.append(helper.make_node(
        "Mul",
        inputs=["token_embeddings", "mask_expanded"],
        outputs=["masked_embeddings"],
    ))

    # 5. Sum across sequence: [batch, hidden_dim]
    nodes.append(helper.make_node(
        "ReduceSum",
        inputs=["masked_embeddings", "reduce_axes"],
        outputs=["summed"],
        keepdims=0,
    ))

    # 6. Sum mask for normalization: [batch, 1]
    nodes.append(helper.make_node(
        "ReduceSum",
        inputs=["mask_expanded", "reduce_axes"],
        outputs=["mask_sum"],
        keepdims=0,
    ))

    # 7. Add epsilon to avoid division by zero
    nodes.append(helper.make_node(
        "Add",
        inputs=["mask_sum", "eps"],
        outputs=["mask_sum_safe"],
    ))

    # 8. Mean pooling: divide by mask sum
    nodes.append(helper.make_node(
        "Div",
        inputs=["summed", "mask_sum_safe"],
        outputs=["pooled"],
    ))

    # 9. Project to final embedding dimension: [batch, embedding_dim]
    nodes.append(helper.make_node(
        "MatMul",
        inputs=["pooled", "proj_weights"],
        outputs=["projected"],
    ))

    # 10. L2 normalize embeddings
    nodes.append(helper.make_node(
        "ReduceL2",
        inputs=["projected"],
        outputs=["norm"],
        axes=[1],
        keepdims=1,
    ))

    nodes.append(helper.make_node(
        "Add",
        inputs=["norm", "eps"],
        outputs=["norm_safe"],
    ))

    nodes.append(helper.make_node(
        "Div",
        inputs=["projected", "norm_safe"],
        outputs=["embeddings"],
    ))

    # Constant tensors for axes
    unsqueeze_axes = numpy_helper.from_array(np.array([2], dtype=np.int64), name="unsqueeze_axes")
    reduce_axes = numpy_helper.from_array(np.array([1], dtype=np.int64), name="reduce_axes")

    # Create the graph
    graph = helper.make_graph(
        nodes=nodes,
        name="test_dna_embedder",
        inputs=[input_ids, attention_mask],
        outputs=[embeddings],
        initializer=[
            token_embed_init,
            proj_init,
            eps_init,
            unsqueeze_axes,
            reduce_axes,
        ],
    )

    # Create the model with opset 13 (compatible with ONNX Runtime 1.20.x)
    model = helper.make_model(
        graph,
        opset_imports=[helper.make_opsetid("", 13)],
        producer_name="llmap-test-generator",
        producer_version="1.0",
        ir_version=8,
    )

    # Add metadata
    model.doc_string = "Test model for LLmap FoundationEmbedder unit tests"

    # Validate
    onnx.checker.check_model(model)

    return model


def main():
    parser = argparse.ArgumentParser(description="Generate test ONNX model for DNA embedding")
    parser.add_argument(
        "--output", "-o",
        default="test_dna_embedder.onnx",
        help="Output path for the ONNX model"
    )
    parser.add_argument(
        "--embedding-dim", "-d",
        type=int,
        default=256,
        help="Output embedding dimension"
    )
    parser.add_argument(
        "--hidden-dim", "-H",
        type=int,
        default=64,
        help="Hidden dimension for token embeddings"
    )
    args = parser.parse_args()

    print(f"Generating test DNA embedder model...")
    print(f"  Embedding dim: {args.embedding_dim}")
    print(f"  Hidden dim: {args.hidden_dim}")

    model = create_test_embedder_model(
        embedding_dim=args.embedding_dim,
        hidden_dim=args.hidden_dim,
    )

    onnx.save(model, args.output)
    print(f"  Saved to: {args.output}")

    # Verify by loading
    loaded = onnx.load(args.output)
    onnx.checker.check_model(loaded)
    print("  Model validation: OK")


if __name__ == "__main__":
    main()
