"""
Fetch FASTA sequences from UniProt REST API.

Uses async HTTP with batching and retries to download protein sequences
for a list of UniProt accession IDs.
"""

import asyncio
import httpx
import time
import os
from typing import List

BATCH_SIZE = 50
MAX_RETRIES = 3
BASE_DELAY = 0.1


def load_entries(filename: str) -> List[str]:
    """Load protein IDs from a text file (one per line)."""
    with open(filename, 'r') as f:
        entries = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(entries)} entry names from {filename}")
    return entries


def _is_valid_fasta(text: str) -> bool:
    return text.strip().startswith('>') and '\n' in text


async def _fetch_one(client: httpx.AsyncClient, entry_name: str, error_log: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{entry_name}.fasta"

    for attempt in range(MAX_RETRIES):
        try:
            resp = await client.get(url)
            if resp.status_code == 200:
                text = resp.text.strip()
                if _is_valid_fasta(text):
                    return text + '\n'
                print(f"Invalid FASTA for {entry_name}")
            elif resp.status_code == 404:
                break
            elif resp.status_code == 429 or resp.status_code >= 500:
                pass  # retry
        except Exception:
            pass

        if attempt < MAX_RETRIES - 1:
            await asyncio.sleep(BASE_DELAY * (2 ** attempt))

    with open(error_log, "a") as f:
        f.write(f"{entry_name}\n")
    return ""


async def _fetch_all(entry_names: List[str], out_fasta: str, error_log: str) -> int:
    with open(out_fasta, "w"):
        pass
    with open(error_log, "w"):
        pass

    async with httpx.AsyncClient(
        timeout=httpx.Timeout(30.0),
        limits=httpx.Limits(max_connections=20, max_keepalive_connections=10),
    ) as client:
        total = len(entry_names)
        for i in range(0, total, BATCH_SIZE):
            batch = entry_names[i:i + BATCH_SIZE]
            tasks = [_fetch_one(client, name, error_log) for name in batch]
            results = await asyncio.gather(*tasks)
            seqs = [r for r in results if r]

            with open(out_fasta, "a") as f:
                f.writelines(seqs)

            print(f"  Batch {i // BATCH_SIZE + 1}: {len(seqs)}/{len(batch)} sequences")

            if i + BATCH_SIZE < total:
                await asyncio.sleep(1.0)

    final_count = sum(1 for line in open(out_fasta) if line.startswith('>'))
    print(f"Total FASTA sequences: {final_count}")
    return final_count


def fetch_fasta(entry_names: List[str], out_fasta: str, error_log: str = "failed_entries.txt") -> int:
    """
    Download FASTA sequences for a list of UniProt IDs.

    Args:
        entry_names: List of UniProt accession IDs.
        out_fasta: Output FASTA file path.
        error_log: Path to log failed IDs.

    Returns:
        Number of sequences successfully downloaded.
    """
    os.makedirs(os.path.dirname(out_fasta) or '.', exist_ok=True)
    t0 = time.perf_counter()
    count = asyncio.run(_fetch_all(entry_names, out_fasta, error_log))
    print(f"Elapsed: {time.perf_counter() - t0:.1f}s")
    return count
